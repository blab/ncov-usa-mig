#!/usr/bin/env Rscript
# Bootstrap CIs for the vaccine period nRR_fixed analysis
# Elderly (65-80y, 81y+) transmission to all age groups, Nov 2020 - May 2021
# Uses same time windows as age_time_RR_analysis.R (2-month rolling, 4-week step)

library(tidyverse)
library(argparse)
library(duckdb)
library(dbplyr)

source("scripts/color_schemes.R")
source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")

normalized_age_rr_fixed <- function(df_rr, baseline_grp, group_vars = NULL){
  if(is.null(group_vars)){
    rr_baseline <- df_rr %>%
      filter(x == baseline_grp & y == baseline_grp) %>%
      select(RR_baseline = RR, N_baseline = N_pairs) %>%
      mutate(RR_baseline = ifelse(N_baseline == 0, NA, RR_baseline))
    df_rr <- df_rr %>%
      mutate(nRR_fixed = RR / rr_baseline$RR_baseline[1])
  } else {
    rr_baseline <- df_rr %>%
      filter(x == baseline_grp & y == baseline_grp) %>%
      select(all_of(group_vars), RR_baseline = RR, N_baseline = N_pairs) %>%
      mutate(RR_baseline = ifelse(N_baseline == 0, NA, RR_baseline))
    df_rr <- df_rr %>%
      left_join(rr_baseline, by = group_vars) %>%
      mutate(nRR_fixed = RR / RR_baseline) %>%
      select(-RR_baseline, -N_baseline)
  }
  return(df_rr)
}

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE,
                      help = "Whether to exclude possible duplicate pairs")
  parser$add_argument('--plot_only', type = 'logical', default = FALSE,
                      help = "Skip bootstrap and reuse existing df_vaccine_period_nRR_ci.tsv files")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
exclude_duplicates <- args$exclude_duplicates
plot_only <- args$plot_only

fn_out_all <- paste0("results/", scenario, "/time_age/df_vaccine_period_nRR_ci.tsv")
fn_out_countries <- paste0("results/", scenario, "/time_age/df_vaccine_period_nRR_ci_countries.tsv")

if (plot_only) {
  message("plot_only mode: reading cached CI TSVs, skipping bootstrap")
  df_vaccine_all <- read_tsv(fn_out_all, show_col_types = FALSE)
  df_vaccine_countries_combined <- read_tsv(fn_out_countries, show_col_types = FALSE)
  countries <- df_vaccine_countries_combined %>% distinct(country) %>% pull(country)
} else {

# Connect to database
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

# Bootstrap parameters
K_BOOT <- 25
SAMP_COV <- 0.8
CI_WIDTH <- 0.95

# Time windows: match age_time_RR_analysis.R parameters
MONTH_BUFFER <- 1
mid_start <- as.Date("2020-03-01")
mid_end <- as.Date("2024-10-01")
all_mid_dates <- seq(from = mid_start, to = mid_end, by = "4 weeks")

# Filter to vaccine period (Nov 2020 - May 2021)
vaccine_mid_dates <- all_mid_dates[all_mid_dates >= as.Date("2020-11-01") &
                                    all_mid_dates <= as.Date("2021-05-31")]
vaccine_lb <- vaccine_mid_dates - 28 * MONTH_BUFFER
vaccine_ub <- vaccine_mid_dates + 28 * MONTH_BUFFER

message(paste0("Processing ", length(vaccine_mid_dates), " time windows for vaccine period"))

# Create country pairs binding
pairs_country <- con %>%
  bind_pairs_exp("country", exclude_duplicates = exclude_duplicates) %>%
  rename(country.x = x, country.y = y) %>%
  mutate(sameCountry = (country.x == country.y))

countries <- tbl(con, "metadata") %>%
  filter(!is.na(country)) %>%
  distinct(country) %>%
  pull(country)

# Point estimates (read from existing TSVs)
fn_series <- paste0("results/", scenario, "/time_age/df_RR_by_time_age_series.tsv")
df_point_all <- read_tsv(fn_series, show_col_types = FALSE) %>%
  filter(date >= as.Date("2020-11-01") & date <= as.Date("2021-05-31"))

fn_countries <- paste0("results/", scenario, "/time_age/df_RR_by_time_age_countries.tsv")
df_point_countries <- read_tsv(fn_countries, show_col_types = FALSE) %>%
  filter(date >= as.Date("2020-11-01") & date <= as.Date("2021-05-31"))

# Bootstrap function for a given set of pairs
run_bootstrap <- function(label, filter_fn = NULL) {
  message(paste0("\n=== Bootstrapping: ", label, " ==="))
  boot_results <- NULL

  for(i in seq_along(vaccine_mid_dates)) {
    message(paste0("Window ", i, "/", length(vaccine_mid_dates), ": ", vaccine_mid_dates[i]))

    boot_rr <- vector("list", K_BOOT)
    for(k in seq_len(K_BOOT)) {
      boot_pairs <- con %>%
        bind_pairs_exp("age_aggr", time_bounds = c(vaccine_lb[i], vaccine_ub[i]),
                       sub_samp = TRUE, samp_cov = SAMP_COV,
                       exclude_duplicates = exclude_duplicates)

      # Apply country filter if provided
      if(!is.null(filter_fn)) {
        boot_pairs <- boot_pairs %>%
          left_join(pairs_country, join_by(strain_1, strain_2)) %>%
          filter_fn()
      }

      boot_iter <- boot_pairs %>%
        calculate_rr_matrix() %>%
        collect() %>%
        normalized_age_rr_fixed(baseline_grp = "26-45y") %>%
        select(x, y, RR, nRR_fixed)

      boot_rr[[k]] <- boot_iter

      if(k %% 10 == 0) {
        gc()
        message(sprintf("  Bootstrap %d/%d", k, K_BOOT))
      }
    }

    boot_combined <- bind_rows(boot_rr)
    ci_vals <- boot_combined %>%
      group_by(x, y) %>%
      summarize(
        RR_lower = quantile(RR, (1 - CI_WIDTH) / 2, na.rm = TRUE),
        RR_upper = quantile(RR, (1 + CI_WIDTH) / 2, na.rm = TRUE),
        nRR_fixed_lower = quantile(nRR_fixed, (1 - CI_WIDTH) / 2, na.rm = TRUE),
        nRR_fixed_upper = quantile(nRR_fixed, (1 + CI_WIDTH) / 2, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(date = vaccine_mid_dates[i])

    if(is.null(boot_results)) {
      boot_results <- ci_vals
    } else {
      boot_results <- bind_rows(boot_results, ci_vals)
    }
  }
  return(boot_results)
}

# All countries combined
boot_all <- run_bootstrap("All Countries")
df_vaccine_all <- df_point_all %>%
  left_join(boot_all, by = c("x", "y", "date"))

write_tsv(df_vaccine_all, fn_out_all)
message(paste0("All countries results saved to ", fn_out_all))

# Country-specific
df_vaccine_countries <- list()
for(ctry in countries) {
  # Create a filter function for this country
  ctry_val <- ctry
  filter_fn <- function(df) df %>% filter(country.x == ctry_val & country.y == ctry_val)

  boot_ctry <- run_bootstrap(ctry, filter_fn)

  df_ctry_point <- df_point_countries %>% filter(country == ctry)
  df_vaccine_countries[[ctry]] <- df_ctry_point %>%
    left_join(boot_ctry, by = c("x", "y", "date"))
}

df_vaccine_countries_combined <- bind_rows(df_vaccine_countries)
write_tsv(df_vaccine_countries_combined, fn_out_countries)
message(paste0("Country results saved to ", fn_out_countries))

DBI::dbDisconnect(con, shutdown = TRUE)
message("Bootstrap complete.")

}  # end if (!plot_only)

# --- Plotting ---
age_order <- c("0-5y", "6-11y", "12-17y", "18-25y", "26-45y", "46-65y", "65-80y", "81y+")

generate_vaccine_ci_plots <- function(df_data, fig_path, title_suffix) {
  plot_data <- df_data %>%
    mutate(
      x = factor(x, levels = age_order),
      y = factor(y, levels = age_order)
    ) %>%
    filter(x %in% c("65-80y", "81y+"))

  # nRR_fixed with CIs
  p_nrr <- ggplot(plot_data, aes(x = date, y = nRR_fixed, color = y, fill = y)) +
    geom_ribbon(aes(ymin = nRR_fixed_lower, ymax = nRR_fixed_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black", alpha = 0.5) +
    facet_wrap(~ x, nrow = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    labs(
      x = "Date",
      y = "nRR",
      color = "To Group",
      fill = "To Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.35, "cm")
    )

  dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
  fn_nrr <- paste0(fig_path, "vaccine_period_nRR_fixed_elderly_ci.png")
  ggsave(fn_nrr, p_nrr, width = 7, height = 3, dpi = 300)
  ggsave(sub("\\.png$", ".svg", fn_nrr), p_nrr, width = 7, height = 3)
  ggsave(sub("\\.png$", "_compact.svg", fn_nrr), p_nrr, width = 5, height = 3)
  message(paste0("nRR plot saved to ", fn_nrr))

  # RR with CIs
  p_rr <- ggplot(plot_data, aes(x = date, y = RR, color = y, fill = y)) +
    geom_ribbon(aes(ymin = RR_lower, ymax = RR_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black", alpha = 0.5) +
    facet_wrap(~ x, nrow = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    scale_y_log10() +
    labs(
      title = paste0("Relative Risk During Vaccine Rollout", title_suffix),
      x = "Date",
      y = "RR (log scale)",
      color = "To Group",
      fill = "To Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )

  fn_rr <- paste0(fig_path, "vaccine_period_RR_elderly_ci.png")
  ggsave(fn_rr, p_rr, width = 14, height = 6, dpi = 300)
  message(paste0("RR plot saved to ", fn_rr))
}

# All countries
fig_path_all <- paste0("figs/", scenario, "/age_time/all_countries/")
generate_vaccine_ci_plots(df_vaccine_all, fig_path_all, " (All Countries)")

# Each country
for(ctry in countries) {
  df_ctry <- df_vaccine_countries_combined %>% filter(country == ctry)
  fig_path_ctry <- paste0("figs/", scenario, "/age_time/", ctry, "/")
  generate_vaccine_ci_plots(df_ctry, fig_path_ctry, paste0(" (", ctry, ")"))
}

cat("\nVaccine period CI analysis complete.\n")
