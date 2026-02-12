library(tidyverse)
library(data.table)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--percentile_threshold', type = 'double', default = 0.98, help = 'Percentile cut off for significant threshold (default is 0.9)')
  parser$add_argument('--min_nb_dist', type = 'integer', default = 2, help = 'Minimum neighbor distance (default is 2 to exclude adjacent states)')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
percentile_threshold <- args$percentile_threshold
min_nb_dist <- args$min_nb_dist

if(percentile_threshold >= 1 | percentile_threshold <= 0){stop("Percentile threhsold should be between 0 and 1!")}
if(min_nb_dist <= 0){stop(("Minimum neighbor distance must be greater than 0!"))}

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

# Filter to only interstate states from the minimum dat
# Exclude pairs where nb_dist < min_nb_dist UNLESS they are international pairs (should always be reported)
inter_state_rr <- state_rr_snap %>%
# filter(region_x != region_y) %>%  # Different regions
  filter(country_x != country_y | is.na(nb_dist) | nb_dist >= min_nb_dist) %>%  # International OR far enough apart
  filter(!is.na(nRR))

nrr_cutoff <- quantile(inter_state_rr$nRR, percentile_threshold, na.rm = TRUE)
message("Cutoff percentile: ",percentile_threshold)
message("nRR Cutoff: ", round(nrr_cutoff,digits = 3))

# Filter to only connections significant from BOTH perspectives
significant_connections <- inter_state_rr %>%
  filter(nRR > nrr_cutoff) %>%
  arrange(date) %>%
  select(date, x, y, nRR, RR, N_pairs, N_x, N_y,
         region_x, region_y, country_x, country_y, nb_dist)

# Save the significant connections
fn_out_path <- paste0("results/", scenario, "/time_state/")
dir.create(file.path(fn_out_path), showWarnings = FALSE)

label_threshold <- round(100*percentile_threshold)
write_tsv(significant_connections,
          file = paste0(fn_out_path, "df_significant_connections_", label_threshold, "_percentile.tsv"))

# Create a summary by time point
summary_by_time <- significant_connections %>%
  group_by(date) %>%
  summarise(
    n_significant_connections = n(),
    nrr_mean = mean(nRR)
  ) %>%
  arrange(date)
