#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)
library(duckdb)
library(dbplyr)

source("scripts/color_schemes.R")
source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")

# Normalization functions (extracted from age_time_RR_analysis.R)
normalized_age_rr_diag <- function(df_rr, group_vars = NULL){
  if(is.null(group_vars)){
    rr_diag <- df_rr %>%
      filter(x == y) %>%
      select(age = x, RR_diag = RR, N_diag = N_pairs) %>%
      mutate(RR_diag = ifelse(N_diag == 0, NA, RR_diag))

    df_rr <- df_rr %>%
      left_join(rr_diag, by = c("x" = "age")) %>%
      rename(RR_xx = RR_diag) %>%
      left_join(rr_diag, by = c("y" = "age")) %>%
      rename(RR_yy = RR_diag) %>%
      mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
  } else {
    rr_diag <- df_rr %>%
      filter(x == y) %>%
      select(all_of(group_vars), age = x, RR_diag = RR, N_diag = N_pairs) %>%
      mutate(RR_diag = ifelse(N_diag == 0, NA, RR_diag))

    df_rr <- df_rr %>%
      left_join(rr_diag, by = c(group_vars, "x" = "age")) %>%
      rename(RR_xx = RR_diag) %>%
      left_join(rr_diag, by = c(group_vars, "y" = "age")) %>%
      rename(RR_yy = RR_diag) %>%
      mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
  }
  return(df_rr)
}

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

# Define time windows
# Pre: August - December 2020
pre_start <- as.Date("2020-08-01")
pre_end <- as.Date("2020-12-31")

# Post: February - June 2021
post_start <- as.Date("2021-02-01")
post_end <- as.Date("2021-06-30")

# Define RR bounds for color scale
UB <- 2.0
LB <- 0.5

# Get metadata and list of states
df_meta <- tbl(con, "metadata")
states <- df_meta %>%
  filter(!is.na(division)) %>%
  distinct(division) %>%
  pull(division)

# Create state pairs binding once
pairs_state <- con %>%
  bind_pairs_exp("division", exclude_duplicates = exclude_duplicates) %>%
  rename(division.x = x, division.y = y) %>%
  mutate(sameDiv = (division.x == division.y))

message("Calculating RR matrices for pre and post periods...")

# Calculate RR matrices for each state and time period
school_state_rr <- NULL

for(time_period in c("Pre", "Post")) {
  if(time_period == "Pre") {
    start_date <- pre_start
    end_date <- pre_end
    period_label <- "Aug-Dec 2020 (Pre)"
  } else {
    start_date <- post_start
    end_date <- post_end
    period_label <- "Feb-Jun 2021 (Post)"
  }

  message(paste0("Processing ", period_label, "..."))

  # Bind pairs for this time period
  pairs_time <- con %>%
    bind_pairs_exp("school_cat", time_bounds = c(start_date, end_date),
                   exclude_duplicates = exclude_duplicates) %>%
    left_join(pairs_state, join_by(strain_1, strain_2))

  for(state in states) {
    # Filter pairs where both strains are from this state
    pair_subset <- pairs_time %>%
      filter(division.x == state & division.y == state) %>%
      collect()

    # Skip if no pairs for this state
    if(nrow(pair_subset) == 0) {
      message(paste0("  Skipping ", state, " (no pairs)"))
      next
    }

    # Calculate RR matrix for this state
    rr_state <- pair_subset %>%
      calculate_rr_matrix() %>%
      mutate(
        state = state,
        time_period = period_label
      ) %>%
      normalized_age_rr_diag() %>%
      normalized_age_rr_fixed(baseline_grp = "Adult")

    # Bind to main dataframe
    if(is.null(school_state_rr)) {
      school_state_rr <- rr_state
    } else {
      school_state_rr <- bind_rows(school_state_rr, rr_state)
    }
  }
}

message("RR matrix calculation complete for individual states.")

# Calculate pooled RR matrices for each trajectory classification
message("Calculating pooled RR matrices for trajectory classifications...")

school_pooled_rr <- NULL

# Get states for each trajectory type
trajectory_states <- df_trajectories %>%
  filter(trajectory_type %in% c("Persistent High", "Persistent Low", "Rising In-Person")) %>%
  select(state, trajectory_type)

for(traj_type in c("Persistent High", "Rising In-Person", "Persistent Low")) {
  message(paste0("Processing pooled ", traj_type, "..."))

  # Get states for this trajectory type
  states_in_traj <- trajectory_states %>%
    filter(trajectory_type == traj_type) %>%
    pull(state)

  for(time_period in c("Pre", "Post")) {
    if(time_period == "Pre") {
      start_date <- pre_start
      end_date <- pre_end
      period_label <- "Aug-Dec 2020 (Pre)"
    } else {
      start_date <- post_start
      end_date <- post_end
      period_label <- "Feb-Jun 2021 (Post)"
    }

    message(paste0("  ", period_label, "..."))

    # Bind pairs for this time period
    pairs_time <- con %>%
      bind_pairs_exp("school_cat", time_bounds = c(start_date, end_date),
                     exclude_duplicates = exclude_duplicates) %>%
      left_join(pairs_state, join_by(strain_1, strain_2))

    # Pool all pairs from states in this trajectory classification
    pair_subset <- pairs_time %>%
      filter(division.x %in% states_in_traj & division.y %in% states_in_traj &
             division.x == division.y) %>%  # Within-state pairs only
      collect()

    # Skip if no pairs
    if(nrow(pair_subset) == 0) {
      message(paste0("  Skipping ", traj_type, " ", time_period, " (no pairs)"))
      next
    }

    # Calculate RR matrix for pooled data
    rr_pooled <- pair_subset %>%
      calculate_rr_matrix() %>%
      mutate(
        trajectory_type = traj_type,
        time_period = period_label,
        pooled = TRUE
      ) %>%
      normalized_age_rr_diag() %>%
      normalized_age_rr_fixed(baseline_grp = "Adult")

    # Bind to main dataframe
    if(is.null(school_pooled_rr)) {
      school_pooled_rr <- rr_pooled
    } else {
      school_pooled_rr <- bind_rows(school_pooled_rr, rr_pooled)
    }
  }
}

# Disconnect from database
DBI::dbDisconnect(con, shutdown = TRUE)

message("Pooled RR matrix calculation complete.")

# Join with trajectory classifications
df_combined <- school_state_rr %>%
  inner_join(df_trajectories %>% select(state, trajectory_type), by = "state") %>%
  filter(trajectory_type %in% c("Persistent High", "Persistent Low", "Rising In-Person"))

# Define age group order
age_order <- c("Pre-School", "Primary School", "Secondary School", "Adult", "Seniors")

df_combined <- df_combined %>%
  mutate(
    x = factor(x, levels = age_order),
    y = factor(y, levels = age_order),
    trajectory_type = factor(trajectory_type,
                            levels = c("Persistent High", "Rising In-Person", "Persistent Low")),
    time_period = factor(time_period,
                        levels = c("Aug-Dec 2020 (Pre)", "Feb-Jun 2021 (Post)"))
  ) %>%
  # Complete the grid to ensure all x-y combinations exist for each state-time_period
  # Group by state, time_period, trajectory_type to avoid creating invalid combinations
  group_by(state, time_period, trajectory_type) %>%
  complete(x, y, fill = list(RR = NA_real_, N_pairs = 0)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(fill_RR = ifelse(N_pairs == 0 | is.na(RR), NA_real_, fill_bound(RR, LB, UB))) %>%
  ungroup()

# Function to create before/after heatmap for a trajectory type
create_trajectory_heatmap <- function(df_data, traj_type, scenario) {

  # Filter to this trajectory type
  df_plot <- df_data %>%
    filter(trajectory_type == traj_type) %>%
    arrange(state, time_period)

  # Count number of states
  n_states <- length(unique(df_plot$state))

  # Create the plot
  p <- ggplot(df_plot, aes(x = x, y = y, fill = fill_RR)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_grid(state ~ time_period, switch = "y") +
    scale_fill_gradient2(
      name = "RR",
      high = RR_HIGH,
      low = RR_LOW,
      limits = c(log10(LB), log10(UB)),
      breaks = c(log10(LB), log10((1+LB)/2), log10(1), log10((1+UB)/2), log10(UB)),
      labels = c(LB, (1+LB)/2, 1, (1+UB)/2, UB),
      na.value = "gray80"
    ) +
    labs(
      title = paste0(traj_type, " States: Before/After RR Heatmaps"),
      subtitle = "Relative Risk (RR) matrices by age group - Academic Year 2020-2021",
      x = "To Age Group",
      y = "From Age Group"
    ) +
    theme_bw() +
    theme(
      strip.text.x = element_text(face = "bold", size = 11),
      strip.text.y = element_text(angle = 180, hjust = 1, size = 9),
      strip.placement = "outside",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank(),
      panel.spacing = unit(0.5, "lines")
    ) +
    coord_equal()

  # Save with height scaled by number of states to maintain same tile size as 3-state plot
  # 3 states at height 8 = ~2.67 inches per state
  # Scale all plots to use the same inches per state
  plot_width <- 10
  plot_height <- n_states * 2.67  # Same vertical space per state for all plots

  fn_out <- paste0("figs/", scenario, "/age_time/before_after_RR_",
                   gsub(" ", "_", tolower(traj_type)), ".png")
  dir.create(dirname(fn_out), recursive = TRUE, showWarnings = FALSE)

  ggsave(fn_out, p, width = plot_width, height = plot_height, dpi = 300)
  message(paste0(traj_type, " heatmap saved to ", fn_out,
                 " (", n_states, " states, height=", plot_height, "in)"))

  return(p)
}

# Create plots for each trajectory type (by state)
trajectory_types <- c("Persistent High", "Rising In-Person", "Persistent Low")

for(traj_type in trajectory_types) {
  create_trajectory_heatmap(df_combined, traj_type, scenario)
}

# Function to create pooled before/after heatmap for a trajectory type
create_pooled_heatmap <- function(df_data, traj_type, scenario) {

  # Define age group order
  age_order <- c("Pre-School", "Primary School", "Secondary School", "Adult", "Seniors")

  # Filter and prepare data
  df_plot <- df_data %>%
    filter(trajectory_type == traj_type) %>%
    mutate(
      x = factor(x, levels = age_order),
      y = factor(y, levels = age_order),
      time_period = factor(time_period,
                          levels = c("Aug-Dec 2020 (Pre)", "Feb-Jun 2021 (Post)"))
    ) %>%
    # Complete the grid to ensure all x-y combinations exist
    group_by(time_period) %>%
    complete(x, y, fill = list(RR = NA_real_, N_pairs = 0)) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(fill_RR = ifelse(N_pairs == 0 | is.na(RR), NA_real_, fill_bound(RR, LB, UB))) %>%
    ungroup()

  # Create the plot
  p <- ggplot(df_plot, aes(x = x, y = y, fill = fill_RR)) +
    geom_tile(color = "white", linewidth = 0.5) +
    facet_grid(. ~ time_period) +
    scale_fill_gradient2(
      name = "RR",
      high = RR_HIGH,
      low = RR_LOW,
      limits = c(log10(LB), log10(UB)),
      breaks = c(log10(LB), log10((1+LB)/2), log10(1), log10((1+UB)/2), log10(UB)),
      labels = c(LB, (1+LB)/2, 1, (1+UB)/2, UB),
      na.value = "gray80"
    ) +
    labs(
      title = paste0(traj_type, " States (Pooled): Before/After RR Heatmaps"),
      subtitle = "Relative Risk (RR) matrices by age group - Academic Year 2020-2021",
      x = "To Age Group",
      y = "From Age Group"
    ) +
    theme_bw() +
    theme(
      strip.text.x = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank(),
      panel.spacing = unit(0.5, "lines")
    ) +
    coord_equal()

  # Save plot
  plot_width <- 8
  plot_height <- 5

  fn_out <- paste0("figs/", scenario, "/age_time/before_after_RR_pooled_",
                   gsub(" ", "_", tolower(traj_type)), ".png")
  dir.create(dirname(fn_out), recursive = TRUE, showWarnings = FALSE)

  ggsave(fn_out, p, width = plot_width, height = plot_height, dpi = 300)
  message(paste0(traj_type, " pooled heatmap saved to ", fn_out))

  return(p)
}

# Create pooled plots for each trajectory type
if(!is.null(school_pooled_rr)) {
  for(traj_type in trajectory_types) {
    create_pooled_heatmap(school_pooled_rr, traj_type, scenario)
  }
}

cat("\nAll before/after heatmaps (by state and pooled) completed.\n")
