#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)
library(duckdb)
library(dbplyr)

source("scripts/color_schemes.R")
source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")

# Normalization: fixed baseline against Adult-Adult RR
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
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
exclude_duplicates <- args$exclude_duplicates

# Connect to database
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

# Read trajectory classifications
fn_trajectories <- paste0("results/", scenario, "/time_age/state_trajectory_classifications.tsv")
df_trajectories <- read_tsv(fn_trajectories, show_col_types = FALSE)

# Define sliding windows: 4-week windows, 1-week step, July 2020 to June 2021
window_size <- 28  # 4 weeks in days
step_size <- 7     # 1 week step

start_range <- as.Date("2020-07-01")
end_range <- as.Date("2021-06-30")

# Generate mid-points
mid_dates <- seq(start_range + window_size / 2, end_range - window_size / 2, by = step_size)
lb_dates <- mid_dates - window_size / 2
ub_dates <- mid_dates + window_size / 2

message(paste0("Processing ", length(mid_dates), " time windows from ",
               min(lb_dates), " to ", max(ub_dates)))

# Create state pairs binding once
pairs_state <- con %>%
  bind_pairs_exp("division", exclude_duplicates = exclude_duplicates) %>%
  rename(division.x = x, division.y = y) %>%
  mutate(sameDiv = (division.x == division.y))

# Get trajectory types and their states
trajectory_types <- c("Persistent High", "Rising In-Person", "Persistent Low", "Other")
trajectory_states <- df_trajectories %>%
  filter(trajectory_type %in% trajectory_types) %>%
  select(state, trajectory_type)

# Bootstrap parameters
K_BOOT <- 25
SAMP_COV <- 0.8
CI_WIDTH <- 0.95

# Calculate pooled nRR_fixed time series for each trajectory type
results <- NULL

for(i in seq_along(mid_dates)) {
  message(paste0("Processing window ", i, "/", length(mid_dates), ": ", mid_dates[i]))

  # Bind pairs for this time window (point estimate)
  pairs_time <- con %>%
    bind_pairs_exp("school_cat", time_bounds = c(lb_dates[i], ub_dates[i]),
                   exclude_duplicates = exclude_duplicates) %>%
    left_join(pairs_state, join_by(strain_1, strain_2))

  for(traj_type in trajectory_types) {
    states_in_traj <- trajectory_states %>%
      filter(trajectory_type == traj_type) %>%
      pull(state)

    # Pool within-state pairs from states in this trajectory
    pair_subset <- pairs_time %>%
      filter(division.x %in% states_in_traj & division.y %in% states_in_traj &
             division.x == division.y) %>%
      collect()

    if(nrow(pair_subset) == 0) {
      message(paste0("  Skipping ", traj_type, " (no pairs)"))
      next
    }

    # Calculate RR matrix and normalize against Adult (point estimate)
    rr_window <- pair_subset %>%
      calculate_rr_matrix() %>%
      mutate(
        date = mid_dates[i],
        trajectory_type = traj_type
      ) %>%
      normalized_age_rr_fixed(baseline_grp = "Adult")

    # Bootstrap CIs
    boot_rr <- vector("list", K_BOOT)
    for(k in seq_len(K_BOOT)) {
      boot_pairs <- con %>%
        bind_pairs_exp("school_cat", time_bounds = c(lb_dates[i], ub_dates[i]),
                       sub_samp = TRUE, samp_cov = SAMP_COV,
                       exclude_duplicates = exclude_duplicates) %>%
        left_join(pairs_state, join_by(strain_1, strain_2)) %>%
        filter(division.x %in% states_in_traj & division.y %in% states_in_traj &
               division.x == division.y) %>%
        collect()

      if(nrow(boot_pairs) == 0) next

      boot_rr[[k]] <- boot_pairs %>%
        calculate_rr_matrix() %>%
        normalized_age_rr_fixed(baseline_grp = "Adult") %>%
        select(x, y, RR, nRR_fixed)
    }

    boot_combined <- bind_rows(boot_rr)
    if(nrow(boot_combined) > 0) {
      ci_vals <- boot_combined %>%
        group_by(x, y) %>%
        summarize(
          RR_lower = quantile(RR, (1 - CI_WIDTH) / 2, na.rm = TRUE),
          RR_upper = quantile(RR, (1 + CI_WIDTH) / 2, na.rm = TRUE),
          nRR_fixed_lower = quantile(nRR_fixed, (1 - CI_WIDTH) / 2, na.rm = TRUE),
          nRR_fixed_upper = quantile(nRR_fixed, (1 + CI_WIDTH) / 2, na.rm = TRUE),
          .groups = "drop"
        )
      rr_window <- rr_window %>%
        left_join(ci_vals, by = c("x", "y"))
    }

    if(is.null(results)) {
      results <- rr_window
    } else {
      results <- bind_rows(results, rr_window)
    }
  }
}

# Disconnect from database
DBI::dbDisconnect(con, shutdown = TRUE)

message("RR calculation complete.")

# Save results
fn_out_tsv <- paste0("results/", scenario, "/time_age/df_nRR_fixed_school_timeseries.tsv")
dir.create(dirname(fn_out_tsv), recursive = TRUE, showWarnings = FALSE)
write_tsv(results, fn_out_tsv)
message(paste0("Results saved to ", fn_out_tsv))

# Filter to interactions involving School Aged
df_plot <- results %>%
  filter(x == "School Aged" | y == "School Aged") %>%
  # Create a label for the interaction partner
  mutate(
    partner = case_when(
      x == "School Aged" & y == "School Aged" ~ "School Aged",
      x == "School Aged" ~ y,
      TRUE ~ x
    )
  ) %>%
  # For off-diagonal, keep only one direction to avoid duplicates
  filter(x == "School Aged") %>%
  mutate(
    partner = factor(partner, levels = c("Pre-School", "School Aged", "Adult", "Seniors")),
    trajectory_type = factor(trajectory_type,
                            levels = c("Persistent High", "Rising In-Person",
                                       "Persistent Low", "Other"))
  )

# nRR_fixed plot
p_nrr <- ggplot(df_plot, aes(x = date, y = nRR_fixed, color = trajectory_type,
                              fill = trajectory_type)) +
  geom_ribbon(aes(ymin = nRR_fixed_lower, ymax = nRR_fixed_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  facet_wrap(~ partner, ncol = 1) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1", name = "Trajectory Type") +
  scale_fill_brewer(palette = "Set1", name = "Trajectory Type") +
  labs(
    title = "School Aged nRR (Adult Baseline) by Interaction Partner",
    subtitle = "4-week sliding windows, July 2020 - June 2021",
    x = "Date",
    y = "nRR (normalized to Adult-Adult, log scale)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

fn_nrr <- paste0("figs/", scenario, "/age_time/school_nRR_fixed_timeseries.png")
dir.create(dirname(fn_nrr), recursive = TRUE, showWarnings = FALSE)
ggsave(fn_nrr, p_nrr, width = 10, height = 12, dpi = 300)
message(paste0("nRR plot saved to ", fn_nrr))

# Raw RR plot
p_rr <- ggplot(df_plot, aes(x = date, y = RR, color = trajectory_type,
                              fill = trajectory_type)) +
  geom_ribbon(aes(ymin = RR_lower, ymax = RR_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  facet_wrap(~ partner, ncol = 1) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1", name = "Trajectory Type") +
  scale_fill_brewer(palette = "Set1", name = "Trajectory Type") +
  labs(
    title = "School Aged RR by Interaction Partner",
    subtitle = "4-week sliding windows, July 2020 - June 2021",
    x = "Date",
    y = "RR (log scale)"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

fn_rr <- paste0("figs/", scenario, "/age_time/school_RR_timeseries.png")
ggsave(fn_rr, p_rr, width = 10, height = 12, dpi = 300)
message(paste0("RR plot saved to ", fn_rr))

cat("\nSchool nRR time series analysis complete.\n")
