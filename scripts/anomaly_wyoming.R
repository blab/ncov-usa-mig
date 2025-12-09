# File: anomaly_wyoming.R
# Author: Amin Bemanian w/ Claude Code
# Date: 2024-11-25
# Description: Investigate Wyoming depleted transmission anomalies with NE states

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(duckdb)
  library(dbplyr)
  library(argparse)
  library(zoo)
  library(tidycensus)
})

source("scripts/color_schemes.R")

#No arguments given this is a specific analysis to a result
scenario <- "CAM_1000"

# Read the time-stratified state RR data
fn_in <- paste0("results/", scenario, "/time_state/df_state_rr_snap.tsv")
state_rr_snap <- fread(fn_in)

# Filter to Wyoming pairs of interest
# Depleted: Maine, Vermont, Nevada, Massachusetts, Rhode Island
# Comparator: Idaho (neighboring state)
state_pairs <- state_rr_snap %>%
  filter(
    (x == "Wyoming" & y %in% c("Maine", "Vermont", "Nevada", "Massachusetts", "Rhode Island", "Idaho")) |
    (y == "Wyoming" & x %in% c("Maine", "Vermont", "Nevada", "Massachusetts", "Rhode Island", "Idaho"))
  ) %>%
  mutate(
    # Create canonical pair labels (Wyoming first)
    pair = case_when(
      x == "Wyoming" ~ paste0("Wyoming-", y),
      y == "Wyoming" ~ paste0("Wyoming-", x),
      TRUE ~ NA_character_
    ),
    # Classify as depleted vs comparator
    pair_type = ifelse(grepl("Idaho", pair), "Comparator (Neighbor)", "Depleted (NE/West)")
  ) %>%
  # Remove duplicates (keep one direction of each pair)
  group_by(date, pair) %>%
  slice_head(n = 1) %>%
  ungroup()

# Create output directory
fn_out_path <- paste0("figs/", scenario, "/anomaly/")
dir.create(fn_out_path, showWarnings = FALSE, recursive = TRUE)

# Define colors for Wyoming pairs
# Use red/orange shades for depleted, green for comparator
pair_colors <- c(
  "Wyoming-Maine" = "#d62728",
  "Wyoming-Vermont" = "#e377c2",
  "Wyoming-Nevada" = "#ff7f0e",
  "Wyoming-Massachusetts" = "#bcbd22",
  "Wyoming-Rhode Island" = "#17becf",
  "Wyoming-Idaho" = "#2ca02c"  # Green for neighbor comparator
)

# Plot RR over time for Wyoming pairs
p_rr_time <- ggplot(state_pairs, aes(x = date, y = RR, color = pair, group = pair)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_y_continuous(
    transform = "log10",
    name = "Relative Risk (RR)",
    breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10),
    labels = c("0.01", "0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10")
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  scale_color_manual(values = pair_colors) +
  labs(
    title = "Relative Risk Over Time: Wyoming Depleted Anomalies",
    subtitle = "Wyoming shows unexpectedly low transmission with NE states; Idaho as neighbor comparator",
    x = "Date",
    color = "State Pair"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save the RR plot
ggsave(paste0(fn_out_path, "wyoming_RR_timeseries.jpg"),
  plot = p_rr_time, width = 12, height = 7, dpi = 192
)

message("RR time series plot saved")

# Plot nRR over time for Wyoming pairs
p_nrr_time <- ggplot(state_pairs, aes(x = date, y = nRR, color = pair, group = pair)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_y_continuous(
    transform = "log10",
    name = "Normalized RR (nRR)",
    breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10),
    labels = c("0.01", "0.05", "0.1", "0.25", "0.5", "1", "2", "5", "10")
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  scale_color_manual(values = pair_colors) +
  labs(
    title = "Normalized RR Over Time: Wyoming Depleted Anomalies",
    x = "Date",
    color = "State Pair"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Save the nRR plot
ggsave(paste0(fn_out_path, "wyoming_nRR_timeseries.jpg"),
  plot = p_nrr_time, width = 12, height = 7, dpi = 192
)

message("nRR time series plot saved")

# ============================================================================
# SEQUENCE COUNTS OVER TIME
# ============================================================================

# Connect to DuckDB for sequence count data
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con_seq <- DBI::dbConnect(duckdb(), fn_db, read_only = TRUE)

# Query sequence counts by state and date
states_of_interest <- c("Wyoming", "Maine", "Vermont", "Nevada", "Massachusetts", "Rhode Island", "Idaho")

df_seq_counts <- tbl(con_seq, "metadata") %>%
  filter(division %in% states_of_interest) %>%
  select(division, date) %>%
  collect() %>%
  mutate(
    date_orig = date,
    date = as.Date(date)
  ) %>%
  # Filter to only full dates (YYYY-MM-DD format)
  filter(nchar(as.character(date_orig)) == 10) %>%
  filter(!is.na(date)) %>%
  filter(date >= as.Date("2020-01-01"), date <= as.Date("2024-12-01")) %>%
  select(division, date)

DBI::dbDisconnect(con_seq, shutdown = TRUE)

# Count sequences by week
df_weekly_counts <- df_seq_counts %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(division, week) %>%
  summarise(n_sequences = n(), .groups = "drop") %>%
  arrange(division, week)

# Create complete grid for all state-week combinations to handle gaps
complete_weeks <- expand.grid(
  division = states_of_interest,
  week = seq(min(df_weekly_counts$week), max(df_weekly_counts$week), by = "week"),
  stringsAsFactors = FALSE
) %>%
  as_tibble()

# Join and fill missing weeks with 0
df_weekly_counts <- complete_weeks %>%
  left_join(df_weekly_counts, by = c("division", "week")) %>%
  mutate(n_sequences = ifelse(is.na(n_sequences), 0, n_sequences)) %>%
  arrange(division, week)

# Apply 4-week simple moving average
df_weekly_counts <- df_weekly_counts %>%
  group_by(division) %>%
  mutate(
    n_sequences_sma = zoo::rollmean(n_sequences, k = 4, fill = NA, align = "center")
  ) %>%
  ungroup()

# Get population data from US Census (2020 Census)
# Use 2020 decennial census for consistency
state_populations <- get_decennial(
  geography = "state",
  variables = "P1_001N", # Total population
  year = 2020,
  survey = "pl"
) %>%
  filter(NAME %in% c("Wyoming", "Maine", "Vermont", "Nevada", "Massachusetts", "Rhode Island", "Idaho")) %>%
  select(NAME, value) %>%
  rename(division = NAME, population = value)

message("\nState populations (2020 Census):")
print(state_populations)

# Join population data and calculate per capita rates
df_weekly_counts <- df_weekly_counts %>%
  left_join(state_populations, by = "division") %>%
  mutate(
    # Sequences per 100,000 population
    n_sequences_per_100k = (n_sequences / population) * 100000,
    n_sequences_sma_per_100k = (n_sequences_sma / population) * 100000
  )

# Define state colors - match pair colors where possible
state_colors <- c(
  "Wyoming" = "#1f77b4",  # Blue for Wyoming (neutral)
  "Maine" = "#d62728",    # Match pair color
  "Vermont" = "#e377c2",  # Match pair color
  "Nevada" = "#ff7f0e",   # Match pair color
  "Massachusetts" = "#bcbd22",  # Match pair color
  "Rhode Island" = "#17becf",   # Match pair color
  "Idaho" = "#2ca02c"     # Match pair color (green)
)

# Plot absolute sequence counts with SMA
p_seq_counts_abs <- ggplot(df_weekly_counts, aes(x = week, y = n_sequences_sma, color = division)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = state_colors, name = "State") +
  scale_y_continuous(
    labels = scales::comma,
    name = "Number of Sequences (4-week SMA)"
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(
    title = "SARS-CoV-2 Weekly Sequence Counts Over Time (Absolute)",
    subtitle = "4-week simple moving average",
    x = "Date"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(paste0(fn_out_path, "sequence_counts_timeseries_absolute.jpg"),
  plot = p_seq_counts_abs, width = 12, height = 7, dpi = 192
)

# Plot per capita sequence counts with SMA
p_seq_counts_per_cap <- ggplot(df_weekly_counts, aes(x = week, y = n_sequences_sma_per_100k, color = division)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = state_colors, name = "State") +
  scale_y_continuous(
    labels = scales::comma,
    name = "Sequences per 100,000 Population (4-week SMA)"
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  labs(
    title = "SARS-CoV-2 Sequencing Effort Over Time (Per Capita)",
    subtitle = "4-week simple moving average normalized by 2020 Census population",
    x = "Date"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(paste0(fn_out_path, "sequence_counts_timeseries_per_capita.jpg"),
  plot = p_seq_counts_per_cap, width = 12, height = 7, dpi = 192
)

message("Sequence counts time series plots saved (absolute and per capita)")

# ============================================================================
# CLADE COMPOSITION OVER TIME
# ============================================================================

# Connect to DuckDB
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db, read_only = TRUE)

# Extract Nextstrain clade information over time for Wyoming and comparison states
states_of_interest <- c("Wyoming", "Maine", "Vermont", "Nevada", "Massachusetts", "Rhode Island", "Idaho")

# Query metadata for these states
df_lineages <- tbl(con, "metadata") %>%
  filter(division %in% states_of_interest) %>%
  select(strain, division, date, Nextstrain_clade) %>%
  collect() %>%
  rename(clade = Nextstrain_clade)

DBI::dbDisconnect(con, shutdown = TRUE)

# Convert date to Date class and filter to reasonable time range
# Only keep sequences with full dates (YYYY-MM-DD format)
df_lineages <- df_lineages %>%
  mutate(
    date_orig = date,
    date = as.Date(date)
  ) %>%
  # Filter to only full dates (length 10: YYYY-MM-DD)
  filter(nchar(as.character(date_orig)) == 10) %>%
  filter(!is.na(date), !is.na(clade), clade != "") %>%
  filter(date >= as.Date("2020-01-01"), date <= as.Date("2024-12-01")) %>%
  mutate(division = factor(division, levels = states_of_interest)) %>%
  select(-date_orig)

# Create monthly bins
df_lineages <- df_lineages %>%
  mutate(
    month = floor_date(date, "month")
  )

# Count sequences by state, month, and clade
lineage_counts <- df_lineages %>%
  group_by(division, month, clade) %>%
  summarise(n = n(), .groups = "drop")

# Get top 20 clades across all data
top_lineages <- lineage_counts %>%
  group_by(clade) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  head(20) %>%
  pull(clade)

message("\nTop 20 Nextstrain clades across all states:")
print(top_lineages)

# Collapse rare clades into "Other" and calculate proportions correctly
lineage_props <- lineage_counts %>%
  # Classify into top clades or "Other"
  mutate(
    lineage_grouped = ifelse(clade %in% top_lineages, clade, "Other")
  ) %>%
  # Group and sum counts for each lineage group
  group_by(division, month, lineage_grouped) %>%
  summarise(n = sum(n), .groups = "drop")

# Create complete grid of all division-month-lineage combinations
# This ensures every lineage appears for every month (even if count = 0)
complete_grid <- expand.grid(
  division = states_of_interest,
  month = unique(lineage_props$month),
  lineage_grouped = c(top_lineages, "Other"),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  mutate(division = factor(division, levels = states_of_interest))

# Join with actual data and fill missing with 0
lineage_props <- complete_grid %>%
  left_join(lineage_props, by = c("division", "month", "lineage_grouped")) %>%
  mutate(n = ifelse(is.na(n), 0, n)) %>%
  # Filter to only keep months where the state actually has data
  group_by(division, month) %>%
  filter(sum(n) > 0) %>%
  # Calculate proportions
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  ungroup()

# Check for any issues with proportions
prop_issues <- lineage_props %>%
  group_by(division, month) %>%
  summarise(total_prop = sum(proportion), .groups = "drop") %>%
  filter(abs(total_prop - 1.0) > 0.01)

if (nrow(prop_issues) > 0) {
  message("WARNING: Some month-state combinations have proportions != 1.0")
  print(head(prop_issues))
}

# Check for any proportion > 1.0
prop_over_1 <- lineage_props %>%
  filter(proportion > 1.0)

if (nrow(prop_over_1) > 0) {
  message("WARNING: Some lineages have proportion > 1.0")
  print(head(prop_over_1, 20))
}


# Create a rainbow color palette with better contrast
# Use a combination of rainbow colors
n_lineages <- length(top_lineages) + 1  # +1 for "Other" (should be 21 total)
lineage_palette <- c(
  scales::hue_pal(h = c(0, 360), l = 55, c = 100)(length(top_lineages)),
  "#888888"  # Gray for "Other"
)
names(lineage_palette) <- c(top_lineages, "Other")

# Area plot for each state
for (state in states_of_interest) {
  df_state <- lineage_props %>%
    filter(division == state) %>%
    # Ensure consistent ordering for stacking
    mutate(lineage_grouped = factor(lineage_grouped, levels = c(top_lineages, "Other"))) %>%
    arrange(month, lineage_grouped) %>%
    # Remove any NA values
    filter(!is.na(proportion), !is.na(month))

  p_lineage <- ggplot(df_state, aes(x = month, y = proportion, fill = lineage_grouped)) +
    geom_area(alpha = 0.8, color = "white", linewidth = 0.2, position = "stack") +
    scale_fill_manual(
      values = lineage_palette,
      name = "Nextstrain Clade"
    ) +
    scale_y_continuous(
      labels = scales::percent,
      expand = c(0, 0)
    ) +
    scale_x_date(
      date_breaks = "6 months",
      date_labels = "%b %Y",
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("Nextstrain Clade Composition Over Time: ", state),
      x = "Date",
      y = "Proportion of Sequences"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  ggsave(paste0(fn_out_path, "clade_composition_", tolower(state), ".jpg"),
    plot = p_lineage, width = 12, height = 7, dpi = 192
  )

  message("Clade composition plot saved for ", state)
}

# Create a combined 7-panel plot (7 rows x 1 column) - Wyoming + 6 comparison states
p_combined <- lineage_props %>%
  # Ensure consistent ordering for stacking
  mutate(lineage_grouped = factor(lineage_grouped, levels = c(top_lineages, "Other"))) %>%
  # Remove any NA values
  filter(!is.na(proportion), !is.na(month)) %>%
  ggplot(aes(x = month, y = proportion, fill = lineage_grouped)) +
  geom_area(alpha = 0.8, color = "white", linewidth = 0.1, position = "stack") +
  scale_fill_manual(
    values = lineage_palette,
    name = "Nextstrain Clade"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = c(0, 0)
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = c(0, 0)
  ) +
  facet_wrap(~division, ncol = 1, scales = "free_x") +
  labs(
    title = "Nextstrain Clade Composition Over Time: Wyoming Anomaly Analysis",
    x = "Date",
    y = "Proportion of Sequences"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  guides(fill = guide_legend(ncol = 5))

ggsave(paste0(fn_out_path, "wyoming_clade_comparison_7panel.jpg"),
  plot = p_combined, width = 14, height = 20, dpi = 192
)

message("Combined 7-panel clade plot saved")

# Summary statistics
message("\n=== Summary Statistics ===")
for (state in states_of_interest) {
  n_seqs <- sum(lineage_counts$n[lineage_counts$division == state])
  date_range <- df_lineages %>%
    filter(division == state) %>%
    summarise(
      min_date = min(date),
      max_date = max(date)
    )
  message(sprintf("%s: %d sequences from %s to %s",
    state, n_seqs, date_range$min_date, date_range$max_date
  ))
}

message("\n=== Analysis Complete ===")
