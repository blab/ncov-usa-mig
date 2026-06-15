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

if (plot_only) {
  message("plot_only mode: reading cached CI TSV, skipping bootstrap")
  df_vaccine_all <- read_tsv(fn_out_all, show_col_types = FALSE)
} else {

# Connect to database
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

# Bootstrap parameters
K_BOOT <- 200
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

# Collect the pair pool once, then bootstrap by subsampling strains in memory
# (same as bind_pairs_exp(sub_samp = TRUE), but without re-querying per resample).
period_lb <- min(vaccine_lb)
period_ub <- max(vaccine_ub)
df_boot_pool <- con %>%
  bind_pairs_exp("age_aggr", time_bounds = c(period_lb, period_ub),
                 exclude_duplicates = exclude_duplicates) %>%
  rename(age.x = x, age.y = y) %>%
  collect()
boot_strains <- unique(c(df_boot_pool$strain_1, df_boot_pool$strain_2))

# Point estimates (read from existing TSV)
fn_series <- paste0("results/", scenario, "/time_age/df_RR_by_time_age_series.tsv")
df_point_all <- read_tsv(fn_series, show_col_types = FALSE) %>%
  filter(date >= as.Date("2020-11-01") & date <= as.Date("2021-05-31"))

# Bootstrap function for a given set of pairs
run_bootstrap <- function(label) {
  message(paste0("\n=== Bootstrapping: ", label, " ==="))
  boot_results <- NULL

  for(i in seq_along(vaccine_mid_dates)) {
    message(paste0("Window ", i, "/", length(vaccine_mid_dates), ": ", vaccine_mid_dates[i]))

    # Pre-filter the in-memory pool to this window
    window_pairs <- df_boot_pool %>%
      filter(transmit_date > vaccine_lb[i] & transmit_date < vaccine_ub[i])

    boot_rr <- map(seq_len(K_BOOT), function(k){
      keep <- sample(boot_strains, floor(length(boot_strains) * SAMP_COV), replace = FALSE)
      window_pairs %>%
        filter(strain_1 %in% keep & strain_2 %in% keep) %>%
        select(strain_1, strain_2, x = age.x, y = age.y) %>%
        calculate_rr_matrix() %>%
        normalized_age_rr_fixed(baseline_grp = "25-64") %>%
        select(x, y, RR, nRR_fixed)
    })

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
set.seed(17)
boot_all <- run_bootstrap("All Countries")
df_vaccine_all <- df_point_all %>%
  left_join(boot_all, by = c("x", "y", "date"))

write_tsv(df_vaccine_all, fn_out_all)
message(paste0("Results saved to ", fn_out_all))

DBI::dbDisconnect(con, shutdown = TRUE)
message("Bootstrap complete.")

}  # end if (!plot_only)

# --- Plotting ---
age_order <- AGE_GROUP_LEVELS  # shared 7-group scheme from color_schemes.R
elderly_groups <- c("65-80", "80+")

generate_vaccine_ci_plots <- function(df_data, fig_path, title_suffix) {
  plot_data <- df_data %>%
    mutate(
      x = factor(x, levels = age_order),
      y = factor(y, levels = age_order)
    ) %>%
    filter(x %in% elderly_groups)

  # nRR_fixed with CIs
  p_nrr <- ggplot(plot_data, aes(x = date, y = nRR_fixed, color = y, fill = y)) +
    geom_ribbon(aes(ymin = nRR_fixed_lower, ymax = nRR_fixed_upper), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black", alpha = 0.5) +
    facet_wrap(~ x, nrow = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
    scale_y_log10() +
    age_group_color_scale(name = "Paired age") +
    age_group_fill_scale(name = "Paired age") +
    labs(
      title = "Elderly RR during Vaccine Rollout",
      x = "Date",
      y = "nRR"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )

  dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
  fn_nrr <- paste0(fig_path, "vaccine_period_nRR_fixed_elderly_ci.png")
  ggsave(fn_nrr, p_nrr, width = 9, height = 4, dpi = 300)
  ggsave(sub("\\.png$", ".svg", fn_nrr), p_nrr, width = 9, height = 4)
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
    age_group_color_scale(name = "Paired age") +
    age_group_fill_scale(name = "Paired age") +
    labs(
      title = paste0("Relative Risk During Vaccine Rollout", title_suffix),
      x = "Date",
      y = "RR (log scale)"
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

fig_path_all <- paste0("figs/", scenario, "/age_time/all_countries/")
generate_vaccine_ci_plots(df_vaccine_all, fig_path_all, "")

cat("\nVaccine period CI analysis complete.\n")
