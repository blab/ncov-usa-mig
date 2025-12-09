library(tidyverse)
library(data.table)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--sd_threshold', type = 'double', default = 2, help = 'Number of standard deviations for significance threshold')
  parser$add_argument('--min_nb_dist', type = 'integer', default = 2, help = 'Minimum neighbor distance (exclude adjacent states)')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
sd_threshold <- args$sd_threshold
min_nb_dist <- args$min_nb_dist

select <- dplyr::select

# Read region data
REGION_DATA <- fread("data/us_states_regions.csv")

# Read neighbor distance data
NB_DIST <- fread("data/nb_dist_states.tsv") %>%
  select(state_x, state_y, nb_dist)

# Read the snapshot RR data (quarterly snapshots)
fn_in <- paste0("results/", scenario, "/time_state/df_state_rr_snap.tsv")
state_rr_snap <- fread(fn_in)

# Add region information
state_rr_snap <- state_rr_snap %>%
  left_join(REGION_DATA %>% select(state, bea_reg, country),
            by = c("x" = "state")) %>%
  rename(region_x = bea_reg, country_x = country) %>%
  left_join(REGION_DATA %>% select(state, bea_reg, country),
            by = c("y" = "state")) %>%
  rename(region_y = bea_reg, country_y = country)

# Add neighbor distance
state_rr_snap <- state_rr_snap %>%
  left_join(NB_DIST, by = c("x" = "state_x", "y" = "state_y"))

# Filter to only inter-state pairs, out-of-region, and non-adjacent states
# Exclude pairs where nb_dist < min_nb_dist UNLESS they are international pairs
# International pairs have country_x != country_y and should always be included
inter_state_rr <- state_rr_snap %>%
  filter(x != y) %>%
  filter(region_x != region_y) %>%  # Different regions
  filter(country_x != country_y | is.na(nb_dist) | nb_dist >= min_nb_dist)  # International OR far enough apart

# For each state (x) and time point (date), calculate mean and SD of out-of-region nRR
out_of_region_stats_x <- inter_state_rr %>%
  group_by(x, date) %>%
  summarise(
    mean_out_region_nRR_x = mean(nRR, na.rm = TRUE),
    sd_out_region_nRR_x = sd(nRR, na.rm = TRUE),
    n_out_region_x = n(),
    .groups = "drop"
  )

# For each state (y) and time point (date), calculate mean and SD of out-of-region nRR
out_of_region_stats_y <- inter_state_rr %>%
  group_by(y, date) %>%
  summarise(
    mean_out_region_nRR_y = mean(nRR, na.rm = TRUE),
    sd_out_region_nRR_y = sd(nRR, na.rm = TRUE),
    n_out_region_y = n(),
    .groups = "drop"
  )

# Join the stats back to the full dataset
inter_state_rr <- inter_state_rr %>%
  left_join(out_of_region_stats_x, by = c("x", "date")) %>%
  left_join(out_of_region_stats_y, by = c("y", "date"))

# Calculate z-score from BOTH states' perspectives
# Z_x = (nRR - mean_x) / SD_x (how unusual is this connection from state x's perspective)
# Z_y = (nRR - mean_y) / SD_y (how unusual is this connection from state y's perspective)
inter_state_rr <- inter_state_rr %>%
  mutate(
    z_score_x = (nRR - mean_out_region_nRR_x) / sd_out_region_nRR_x,
    z_score_y = (nRR - mean_out_region_nRR_y) / sd_out_region_nRR_y,
    is_significant_x = abs(z_score_x) > sd_threshold,
    is_significant_y = abs(z_score_y) > sd_threshold,
    is_significant_both = is_significant_x & is_significant_y
  )

# Filter to only connections significant from BOTH perspectives
significant_connections <- inter_state_rr %>%
  filter(is_significant_both) %>%
  arrange(date, desc(pmin(abs(z_score_x), abs(z_score_y)))) %>%
  select(date, x, y, nRR, RR, N_pairs, N_x, N_y,
         region_x, region_y, country_x, country_y, nb_dist,
         mean_out_region_nRR_x, sd_out_region_nRR_x, z_score_x, n_out_region_x,
         mean_out_region_nRR_y, sd_out_region_nRR_y, z_score_y, n_out_region_y)

# Save the significant connections
fn_out_path <- paste0("results/", scenario, "/time_state/")
dir.create(file.path(fn_out_path), showWarnings = FALSE)

write_tsv(significant_connections,
          file = paste0(fn_out_path, "df_significant_connections_", sd_threshold, "sd.tsv"))

# Create a summary by time point
summary_by_time <- significant_connections %>%
  group_by(date) %>%
  summarise(
    n_significant_connections = n(),
    mean_z_score_x = mean(abs(z_score_x)),
    mean_z_score_y = mean(abs(z_score_y)),
    max_z_score_x = max(abs(z_score_x)),
    max_z_score_y = max(abs(z_score_y)),
    states_involved = n_distinct(c(x, y)),
    .groups = "drop"
  ) %>%
  arrange(date)

write_tsv(summary_by_time,
          file = paste0(fn_out_path, "summary_significant_connections_by_time_", sd_threshold, "sd.tsv"))

# Print summary
cat("\n=== Significant Connection Analysis Summary ===\n")
cat("Threshold:", sd_threshold, "standard deviations (from BOTH states' perspectives)\n")
cat("Minimum neighbor distance:", min_nb_dist, "(adjacent states excluded)\n")
cat("Data source: Quarterly snapshots (3-month windows)\n")
cat("Total significant connections found:", nrow(significant_connections), "\n")
cat("Time points with significant connections:", n_distinct(significant_connections$date), "\n")
cat("\nTop 10 most significant connections:\n")
print(head(significant_connections %>%
           select(date, x, y, nRR, nb_dist, z_score_x, z_score_y) %>%
           arrange(desc(pmin(abs(z_score_x), abs(z_score_y)))), 10))

cat("\n=== Summary by time point ===\n")
print(summary_by_time)
