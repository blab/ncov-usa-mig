# File: ohio_oklahoma_anomaly.R
# Author: Amin Bemanian w/ Claude Code
# Date: 2024-11-25
# Description: Investigate Ohio-Oklahoma transmission anomaly

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(duckdb)
  library(dbplyr)
  library(argparse)
  library(zoo)
  library(tidycensus)
  library(igraph)
})

source("scripts/color_schemes.R")
source("scripts/bind_pairs_exp.R")

#No arguments given this is a specific analysis to a result
scenario <- "CAM_1000"

# Read the time-stratified state RR data
fn_in <- paste0("results/", scenario, "/time_state/df_state_rr_snap.tsv")
state_rr_snap <- fread(fn_in)

# Filter to the three state pairs of interest
state_pairs <- state_rr_snap %>%
  filter(
    (x == "Ohio" & y == "Oklahoma") |
    (x == "Oklahoma" & y == "Ohio") |
    (x == "Ohio" & y == "Illinois") |
    (x == "Illinois" & y == "Ohio") |
    (x == "Oklahoma" & y == "Texas") |
    (x == "Texas" & y == "Oklahoma")
  ) %>%
  mutate(
    # Create canonical pair labels (alphabetically ordered)
    pair = case_when(
      (x == "Ohio" & y == "Oklahoma") | (x == "Oklahoma" & y == "Ohio") ~ "Ohio-Oklahoma",
      (x == "Ohio" & y == "Illinois") | (x == "Illinois" & y == "Ohio") ~ "Illinois-Ohio",
      (x == "Oklahoma" & y == "Texas") | (x == "Texas" & y == "Oklahoma") ~ "Oklahoma-Texas",
      TRUE ~ NA_character_
    )
  ) %>%
  # Remove duplicates (keep one direction of each pair)
  group_by(date, pair) %>%
  slice_head(n = 1) %>%
  ungroup()

# Create output directory
fn_out_path <- paste0("figs/", scenario, "/anomaly/")
dir.create(fn_out_path, showWarnings = FALSE, recursive = TRUE)

# Define consistent colors for state pairs
pair_colors <- c(
  "Ohio-Oklahoma" = "firebrick",
  "Illinois-Ohio" = "steelblue",
  "Oklahoma-Texas" = "forestgreen"
)

# Plot RR over time for the three pairs
p_rr_time <- ggplot(state_pairs, aes(x = date, y = RR, color = pair, group = pair)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_y_continuous(
    transform = "log10",
    name = "Relative Risk (RR)",
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
    labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10")
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  scale_color_manual(values = pair_colors) +
  labs(
    title = "Relative Risk Over Time: Ohio-Oklahoma vs Neighboring States",
    subtitle = "Raw RR values before normalization",
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
ggsave(paste0(fn_out_path, "ohio_oklahoma_RR_timeseries.jpg"),
  plot = p_rr_time, width = 12, height = 7, dpi = 192
)

message("RR time series plot saved to: ", fn_out_path, "ohio_oklahoma_RR_timeseries.jpg")

# Plot nRR over time for the three pairs
p_nrr_time <- ggplot(state_pairs, aes(x = date, y = nRR, color = pair, group = pair)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  scale_y_continuous(
    transform = "log10",
    name = "Normalized RR (nRR)",
    breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10),
    labels = c("0.1", "0.25", "0.5", "1", "2", "5", "10")
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b %Y") +
  scale_color_manual(values = pair_colors) +
  labs(
    title = "Normalized RR Over Time: Ohio-Oklahoma vs Neighboring States",
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
ggsave(paste0(fn_out_path, "ohio_oklahoma_nRR_timeseries.jpg"),
  plot = p_nrr_time, width = 12, height = 7, dpi = 192
)

message("nRR time series plot saved to: ", fn_out_path, "ohio_oklahoma_nRR_timeseries.jpg")

# ============================================================================
# SEQUENCE COUNTS OVER TIME
# ============================================================================

# Connect to DuckDB for sequence count data
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con_seq <- DBI::dbConnect(duckdb(), fn_db, read_only = TRUE)

# Query sequence counts by state and date
states_of_interest <- c("Ohio", "Oklahoma", "Illinois", "Texas")

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
  filter(NAME %in% c("Ohio", "Oklahoma", "Illinois", "Texas")) %>%
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

# Define state colors
state_colors <- c(
  "Ohio" = "firebrick",
  "Oklahoma" = "darkorange",
  "Illinois" = "steelblue",
  "Texas" = "forestgreen"
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

# Extract Nextstrain clade information over time for the four states
states_of_interest <- c("Ohio", "Oklahoma", "Illinois", "Texas")

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

# Create a combined 4-panel plot (4 rows x 1 column)
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
    title = "Nextstrain Clade Composition Over Time",
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

ggsave(paste0(fn_out_path, "ohio_oklahoma_clade_comparison_4panel.jpg"),
  plot = p_combined, width = 14, height = 16, dpi = 192
)

message("Combined 4-panel clade plot saved")

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

# ============================================================================
# NETWORK OF RELATED PAIRS
# ============================================================================

#Obtain strains for the pairs of Ohio-Oklahoma identical sequences
oh_ok_strain_names <- bind_pairs_exp(con_seq, "division") %>%
  filter((x == "Ohio" & y == "Oklahoma") |
    (x == "Oklahoma" & y == "Ohio")
  ) %>% select(strain_1,strain_2) %>%
  {
    union_all(
      select(., strain = strain_1),
      select(., strain = strain_2)
    )
  } %>%
  distinct(strain) %>%
  pull(strain)

#Isolate out the OH-OK specific sequences to compare against the general Ohio and Oklahoma sequences
df_oh_ok_seq <- tbl(con_seq,"metadata") %>%
  filter(strain %in% oh_ok_strain_names) %>%
  collect()

#Find other sequences that are identical outside of the OH-OK ones
connected_pair_list <- bind_pairs_exp(con_seq,"division") %>%
  filter(strain_1 %in% oh_ok_strain_names | strain_2 %in% oh_ok_strain_names) %>%
  select(strain_1,strain_2) %>%
  collect()

connected_pair_strain_names <- bind_pairs_exp(con_seq,"division") %>%
  filter(strain_1 %in% oh_ok_strain_names |strain_2 %in% oh_ok_strain_names) %>%
  select(strain_1,strain_2) %>%
  {
    union_all(
      select(., strain = strain_1),
      select(., strain = strain_2)
    )
  } %>%
  distinct(strain) %>%
  pull(strain)

df_connected_seq <- tbl(con_seq,"metadata") %>%
  filter(strain %in% connected_pair_strain_names) %>%
  collect()

message("Counts of sequences identical to ones involved in Ohio-Oklahoma pairs")
message(sprintf("Number of sequences in OH or OK: %d \nNumber of sequences across all states: %d",
  df_oh_ok_seq %>% nrow(),
  df_connected_seq %>% nrow
))

DBI::dbDisconnect(con_seq, shutdown = TRUE)

# ============================================================================
# NETWORK ANALYSIS: BUILD GRAPH AND FIND MAXIMAL CONENCTED COMPONENTS
# ============================================================================

message("\n=== Building network graph from connected pairs ===")

# Create igraph object from edge list
g <- graph_from_data_frame(
  d = connected_pair_list,
  directed = FALSE,
  vertices = data.frame(name = connected_pair_strain_names)
)

message(sprintf("Graph created: %d nodes, %d edges", vcount(g), ecount(g)))


components <- components(g)
message(sprintf("Found %d connected components", components$no))

# Add component ID to df_connected_seq
component_membership <- data.frame(
  strain = names(components$membership),
  component_id = as.integer(components$membership),
  row.names = NULL
)

df_connected_seq <- df_connected_seq %>%
  left_join(component_membership, by = "strain")

# Calculate first and last date for each component
component_dates <- df_connected_seq %>%
  mutate(
    date_orig = date,
    date = as.Date(date)
  ) %>%
  # Filter to only full dates (YYYY-MM-DD format, length 10)
  filter(nchar(as.character(date_orig)) == 10) %>%
  filter(!is.na(date), !is.na(component_id)) %>%
  group_by(component_id) %>%
  summarise(
    first_date = min(date),
    last_date = max(date),
    date_span_days = as.numeric(max(date) - min(date)),
    .groups = "drop"
  )

# Join date info back to component_membership
component_membership <- component_membership %>%
  left_join(component_dates, by = "component_id")

# Summary of component sizes
message(sprintf("Largest component has %d nodes", max(components$csize)))
message(sprintf("Smallest component has %d nodes", min(components$csize)))
message(sprintf("Number of singleton components (size=1): %d", sum(components$csize == 1)))

# Investigate singletons if they exist
if (sum(components$csize == 1) > 0) {
  message("\n=== INVESTIGATING SINGLETON COMPONENTS ===")
  singleton_ids <- which(components$csize == 1)
  singleton_strains <- names(components$membership)[components$membership %in% singleton_ids]

  message(sprintf("Found %d singleton components with %d total sequences",
                  length(singleton_ids), length(singleton_strains)))

  # Check if these strains appear in the edge list
  singletons_in_edges <- connected_pair_list %>%
    filter(strain_1 %in% singleton_strains | strain_2 %in% singleton_strains)

  message(sprintf("Singleton strains appearing in edge list: %d edges", nrow(singletons_in_edges)))

  if (nrow(singletons_in_edges) > 0) {
    message("WARNING: Singleton strains found in edge list - this indicates a graph construction issue!")
    message("First few singleton edges:")
    print(head(singletons_in_edges, 10))
  } else {
    message("Singleton strains do NOT appear in edge list")
    message("This suggests isolated vertices were added during graph construction")
  }

  # Check metadata for singleton sequences
  singleton_metadata <- df_connected_seq %>%
    filter(strain %in% singleton_strains) %>%
    select(strain, division, date, Nextstrain_clade) %>%
    head(20)

  message("\nMetadata for first 20 singleton sequences:")
  print(singleton_metadata)
}

# Create histogram of component sizes
# Pre-bin the data to avoid empty bins
component_sizes_df <- data.frame(size = components$csize)

# Create bins manually and count
breaks <- seq(0, max(component_sizes_df$size), length.out = 51)
component_hist <- hist(component_sizes_df$size, breaks = breaks, plot = FALSE)

# Create data frame from histogram, filtering out zero counts
hist_df <- data.frame(
  midpoint = component_hist$mids,
  count = component_hist$counts
) %>%
  filter(count > 0)  # Remove bins with zero counts

p_component_sizes <- ggplot(hist_df, aes(x = midpoint, y = count)) +
  geom_col(fill = "steelblue", color = "white", width = diff(breaks)[1] * 0.9) +
  scale_y_continuous(
    name = "Number of Components",
    limits=c(0,NA)
  ) +
  scale_x_continuous(
    name = "Component Size (number of sequences)",
    labels = scales::comma
  ) +
  labs(
    title = "Distribution of Connected Component Sizes",
    subtitle = sprintf("%d components total, largest has %d nodes",
                       components$no, max(components$csize))
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(paste0(fn_out_path, "component_size_distribution.jpg"),
  plot = p_component_sizes, width = 10, height = 7, dpi = 192
)

message("Component size histogram saved")

# ============================================================================
# SANITY CHECK: Composition of size-2 components
# ============================================================================

message("\n=== Analyzing size-2 components ===")

# Get all size-2 component IDs
size2_comp_ids <- which(components$csize == 2)
message(sprintf("Number of size-2 components: %d", length(size2_comp_ids)))

# For each size-2 component, get the two states involved
size2_compositions <- lapply(size2_comp_ids, function(comp_id) {
  strains_in_comp <- names(components$membership)[components$membership == comp_id]

  states <- df_connected_seq %>%
    filter(strain %in% strains_in_comp) %>%
    pull(division) %>%
    sort()

  if (length(states) == 2) {
    return(data.frame(
      component_id = comp_id,
      state_1 = states[1],
      state_2 = states[2],
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
}) %>%
  bind_rows()

# Create state-pair labels
size2_compositions <- size2_compositions %>%
  mutate(
    state_pair = paste(state_1, state_2, sep = "-"),
    is_oh_ok = (state_1 == "Ohio" & state_2 == "Oklahoma") |
               (state_1 == "Oklahoma" & state_2 == "Ohio")
  )

# Count by state pair type
size2_summary <- size2_compositions %>%
  group_by(state_pair) %>%
  summarise(n_components = n(), .groups = "drop") %>%
  arrange(desc(n_components))

message("\nSize-2 component composition by state pairs:")
print(size2_summary)

# Count OH-OK vs others
n_oh_ok_pairs <- sum(size2_compositions$is_oh_ok)
n_other_pairs <- nrow(size2_compositions) - n_oh_ok_pairs
pct_oh_ok <- (n_oh_ok_pairs / nrow(size2_compositions)) * 100

message(sprintf("\nOH-OK pairs: %d (%.1f%%)", n_oh_ok_pairs, pct_oh_ok))
message(sprintf("Other state pairs: %d (%.1f%%)", n_other_pairs, 100 - pct_oh_ok))

# Show some examples of non-OH-OK size-2 components
if (n_other_pairs > 0) {
  message("\nExamples of non-OH-OK size-2 components:")
  non_oh_ok_examples <- size2_compositions %>%
    filter(!is_oh_ok) %>%
    head(10)
  print(non_oh_ok_examples)
}

message("Component size histogram saved")

# Plot components over time (first date vs size)
component_time_data <- data.frame(
  component_id = seq_along(components$csize),
  size = components$csize
) %>%
  left_join(component_dates, by = "component_id") %>%
  filter(!is.na(first_date))

p_component_time <- ggplot(component_time_data, aes(x = first_date, y = size)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 2) +
  scale_y_continuous(
    name = "Component Size (number of sequences)",
    labels = scales::comma
  ) +
  scale_x_date(
    name = "First Sequence Date",
    date_breaks = "6 months",
    date_labels = "%b %Y",
    limits = c(as.Date("2020-01-01"),as.Date("2024-10-31"))
  ) +
  labs(
    title = "Connected Components Over Time",
    subtitle = sprintf("%d components with temporal data", nrow(component_time_data))
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(fn_out_path, "component_size_over_time.jpg"),
  plot = p_component_time, width = 12, height = 7, dpi = 192
)

message("Component size over time plot saved")

# ============================================================================
# MAP OF LARGEST CONNECTED COMPONENT
# ============================================================================

source("scripts/cam_map.R")

# Identify the largest component
largest_comp_id <- which.max(components$csize)
message(sprintf("\n=== Mapping largest component (ID: %d, size: %d) ===",
                largest_comp_id, components$csize[largest_comp_id]))

# Get sequences in the largest component
df_largest_comp <- df_connected_seq %>%
  filter(component_id == largest_comp_id) %>%
  filter(!is.na(division))

# Count sequences per state in the largest component
state_counts <- df_largest_comp %>%
  group_by(division) %>%
  summarise(n_sequences = n(), .groups = "drop")

message(sprintf("Largest component spans %d states", nrow(state_counts)))
print(state_counts)

# Get edges within the largest component, aggregated by state pairs
strains_in_comp <- df_largest_comp$strain

edges_in_comp <- connected_pair_list %>%
  filter(strain_1 %in% strains_in_comp & strain_2 %in% strains_in_comp)

# Join with state information
edges_with_states <- edges_in_comp %>%
  left_join(df_largest_comp %>% select(strain, division_1 = division),
            by = c("strain_1" = "strain")) %>%
  left_join(df_largest_comp %>% select(strain, division_2 = division),
            by = c("strain_2" = "strain")) %>%
  filter(!is.na(division_1), !is.na(division_2))

# Aggregate edges by state pairs (undirected)
state_edges <- edges_with_states %>%
  mutate(
    state_x = pmin(division_1, division_2),
    state_y = pmax(division_1, division_2)
  ) %>%
  group_by(state_x, state_y) %>%
  summarise(n_edges = n(), .groups = "drop") %>%
  filter(state_x != state_y)  # Remove self-loops

message(sprintf("Found %d inter-state connections", nrow(state_edges)))

# Load CAM map
cam_map <- prep_cam_map()

# Get state centroids for mapping
state_centroids <- cam_map %>%
  filter(NAME_En != "Hawaii") %>%
  st_point_on_surface() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(division = cam_map %>% filter(NAME_En != "Hawaii") %>% pull(NAME_En)) %>%
  rename(longitude = X, latitude = Y)

# Join counts with centroids
state_map_data <- state_counts %>%
  inner_join(state_centroids, by = "division")

# Get date range for the largest component
largest_comp_dates <- component_dates %>%
  filter(component_id == largest_comp_id)

# Create the network map
p_network_map <- ggplot() +
  # Base CAM map
  geom_sf(data = cam_map %>% filter(NAME_En != "Hawaii"),
          fill = "gray95", color = "gray70", linewidth = 0.3) +
  # State nodes (circles sized by sequence count)
  geom_point(data = state_map_data,
             aes(x = longitude, y = latitude, size = n_sequences),
             color = "steelblue", alpha = 0.7) +
  scale_size_area(name = "Sequences", max_size = 15) +
  scale_alpha_continuous(name = "Connections", range = c(0.2, 0.8)) +
  labs(
    title = "Geography of Largest Connected Component",
    subtitle = sprintf("%d sequences across %d states/provinces\nDate range: %s to %s",
                      sum(state_counts$n_sequences), nrow(state_counts),
                      largest_comp_dates$first_date, largest_comp_dates$last_date)
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(paste0(fn_out_path, "largest_component_network_map.jpg"),
  plot = p_network_map, width = 14, height = 10, dpi = 192
)

message("Network map saved")

# ============================================================================
# STATE INVOLVEMENT IN CONNECTED COMPONENTS
# ============================================================================

message("\n=== Analyzing state involvement in connected components ===")

# For each component, get the list of states involved
component_states <- df_connected_seq %>%
  filter(!is.na(component_id), !is.na(division)) %>%
  group_by(component_id) %>%
  summarise(
    states = list(unique(division)),
    n_states = n_distinct(division),
    n_sequences = n(),
    .groups = "drop"
  )

# Total number of components
total_components <- length(unique(component_states$component_id))
message(sprintf("Total components with state data: %d", total_components))

# For each state (excluding OH and OK), count how many components it appears in
all_states <- df_connected_seq %>%
  filter(!is.na(division)) %>%
  pull(division) %>%
  unique() %>%
  setdiff(c("Ohio", "Oklahoma"))

state_involvement <- sapply(all_states, function(state) {
  sum(sapply(component_states$states, function(states_in_comp) {
    state %in% states_in_comp
  }))
})

state_involvement_df <- data.frame(
  state = all_states,
  n_components = state_involvement,
  pct_components = (state_involvement / total_components) * 100
) %>%
  arrange(desc(pct_components))

message("\nTop 10 states by component involvement:")
print(head(state_involvement_df, 10))

# Calculate percentage of components that ONLY have OH and OK (no outside states)
oh_ok_only_components <- component_states %>%
  filter(sapply(states, function(s) {
    all(s %in% c("Ohio", "Oklahoma"))
  }))

pct_oh_ok_only <- (nrow(oh_ok_only_components) / total_components) * 100
message(sprintf("\nComponents with ONLY Ohio and Oklahoma: %d (%.1f%%)",
                nrow(oh_ok_only_components), pct_oh_ok_only))

# Prepare data for barplot: top 10 states + "OH+OK only" category
plot_data <- state_involvement_df %>%
  head(10) %>%
  bind_rows(
    data.frame(
      state = "OH + OK only",
      n_components = nrow(oh_ok_only_components),
      pct_components = pct_oh_ok_only
    )
  ) %>%
  # Create factor with levels in reverse order for top-to-bottom plotting
  mutate(state = factor(state, levels = rev(state)))

# Create vertical barplot
p_state_involvement <- ggplot(plot_data, aes(x = state, y = pct_components)) +
  geom_col(
    aes(fill = state == "OH + OK only"),
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_components)),
    hjust = -0.2,
    size = 3.5
  ) +
  scale_fill_manual(values = c("steelblue", "darkorange")) +
  scale_y_continuous(
    name = "Percentage of Components (%)",
    limits = c(0, max(plot_data$pct_components) * 1.15),
    expand = c(0, 0)
  ) +
  labs(
    title = "State Involvement in OH-OK Connected Components",
    subtitle = sprintf("Analysis of %d connected components", total_components),
    x = NULL
  ) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold")
  )

ggsave(paste0(fn_out_path, "state_involvement_barplot.jpg"),
  plot = p_state_involvement, width = 10, height = 8, dpi = 192
)

message("\nState involvement barplot saved to: ", fn_out_path, "state_involvement_barplot.jpg")

# Save the state involvement data
fn_state_involvement <- paste0("results/", scenario, "/anomaly/state_involvement_in_components.tsv")
dir.create(dirname(fn_state_involvement), showWarnings = FALSE, recursive = TRUE)
fwrite(state_involvement_df, fn_state_involvement, sep = "\t")

message("State involvement data saved to: ", fn_state_involvement)

message("\n=== Network Analysis Complete ===")