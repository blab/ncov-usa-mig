#File: age_season_fl_mideast.R
#Author(s): Amin Bemanian
#Date: 6/9/26
#Description: Age x age RR matrices for identical-sequence pairs restricted to
# Florida and the Mideast (BEA region), stratified by geographic pairing
# (FL/FL, Mideast/Mideast, FL/Mideast) and by meteorological season
# (Mar-May, Jun-Aug, Sep-Nov, Dec-Feb), pooled across all years.
# Produces a 4 (season) x 3 (geography) grid of RR heatmaps.

library(argparse)
library(tidyverse)
library(duckdb)
library(dbplyr)
library(lubridate)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/color_schemes.R")

# RR heatmap color bounds (linear RR terms)
LB <- 0.5
UB <- 2.0

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--bin_size', type = 'integer', default = 5, help = 'Age bin width in years (default 5)')
  parser$add_argument('--age_cap', type = 'integer', default = 85, help = 'Top-coded age; ages >= cap pooled into a single bin (default 85)')
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE, help = "Whether to exclude possible duplicate pairs, default is FALSE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
bin_size <- args$bin_size
age_cap <- args$age_cap
exclude_duplicates <- args$exclude_duplicates

# Configurable, top-coded age binning
age_bin <- function(a, bin = 5, cap = 85){
  ifelse(a >= cap,
         sprintf("%d+y", cap),
         sprintf("%02d-%02dy", (a %/% bin) * bin, (a %/% bin) * bin + (bin - 1)))
}

# Ordered factor levels for the age bins (drives heatmap axis order)
age_bin_levels <- function(bin = 5, cap = 85){
  starts <- seq(0, cap - bin, by = bin)
  c(sprintf("%02d-%02dy", starts, starts + (bin - 1)), sprintf("%d+y", cap))
}

SEASON_LEVELS <- c("Mar-May", "Jun-Aug", "Sep-Nov", "Dec-Feb")
GEO_LEVELS <- c("FL/FL", "Mideast/Mideast", "FL/Mideast")

assign_season <- function(m){
  case_when(
    m %in% c(3, 4, 5)   ~ "Mar-May",
    m %in% c(6, 7, 8)   ~ "Jun-Aug",
    m %in% c(9, 10, 11) ~ "Sep-Nov",
    TRUE                ~ "Dec-Feb"
  )
}

fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

# Full date range so the pairs_time table is used and transmit_date is carried
full_bounds <- c(as.Date("2019-01-01"), as.Date("2025-01-01"))

# Bind division and bea_reg to classify geography. pairs_div carries
# transmit_date through from the pairs_time table. Age is joined separately
# below because bind_pairs_exp's character "NA" filter is incompatible with the
# integer age_adj column.
pairs_div <- con %>%
  bind_pairs_exp("division", time_bounds = full_bounds, exclude_duplicates = exclude_duplicates) %>%
  rename(division.x = x, division.y = y)

pairs_reg <- con %>%
  bind_pairs_exp("bea_reg", time_bounds = full_bounds, exclude_duplicates = exclude_duplicates) %>%
  rename(bea_reg.x = x, bea_reg.y = y)

# Local strain -> raw age dictionary
age_dict <- tbl(con, "metadata") %>%
  select(strain, age_adj) %>%
  collect()

df_pairs <- pairs_div %>%
  left_join(pairs_reg %>% select(strain_1, strain_2, bea_reg.x, bea_reg.y), by = join_by(strain_1, strain_2)) %>%
  collect() %>%
  left_join(age_dict %>% rename(age.x = age_adj), by = join_by(strain_1 == strain)) %>%
  left_join(age_dict %>% rename(age.y = age_adj), by = join_by(strain_2 == strain)) %>%
  filter(!is.na(age.x) & !is.na(age.y))

# Classify each strain as FL / Mideast / NA, restrict to FL+Mideast pairs
df_pairs <- df_pairs %>%
  mutate(
    geo.x = case_when(division.x == "Florida" ~ "FL", bea_reg.x == "Mideast" ~ "Mideast", TRUE ~ NA_character_),
    geo.y = case_when(division.y == "Florida" ~ "FL", bea_reg.y == "Mideast" ~ "Mideast", TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(geo.x) & !is.na(geo.y)) %>%
  mutate(
    geo_cat = case_when(
      geo.x == "FL" & geo.y == "FL" ~ "FL/FL",
      geo.x == "Mideast" & geo.y == "Mideast" ~ "Mideast/Mideast",
      TRUE ~ "FL/Mideast"
    ),
    geo_cat = factor(geo_cat, levels = GEO_LEVELS),
    season = factor(assign_season(month(transmit_date)), levels = SEASON_LEVELS)
  )

# Apply age binning; x/y become the exposure categories for calculate_rr_matrix
bin_levels <- age_bin_levels(bin_size, age_cap)
df_pairs <- df_pairs %>%
  mutate(
    x = factor(age_bin(age.x, bin_size, age_cap), levels = bin_levels),
    y = factor(age_bin(age.y, bin_size, age_cap), levels = bin_levels)
  ) %>%
  filter(!is.na(x) & !is.na(y))

# Compute RR matrix per (season x geo_cat) stratum
strata <- df_pairs %>% distinct(season, geo_cat)

df_rr <- pmap_dfr(strata, function(season, geo_cat){
  subset <- df_pairs %>%
    filter(season == !!season & geo_cat == !!geo_cat) %>%
    mutate(x = as.character(x), y = as.character(y)) %>%
    select(strain_1, strain_2, x, y)
  if(nrow(subset) == 0) return(NULL)
  subset %>%
    calculate_rr_matrix() %>%
    mutate(season = season, geo_cat = geo_cat)
})

df_rr <- df_rr %>%
  mutate(
    season = factor(season, levels = SEASON_LEVELS),
    geo_cat = factor(geo_cat, levels = GEO_LEVELS),
    x = factor(x, levels = bin_levels),
    y = factor(y, levels = bin_levels)
  )

# Write table
fn_out_path <- paste0("results/", scenario, "/season_geo_age/")
dir.create(fn_out_path, recursive = TRUE, showWarnings = FALSE)
write_tsv(df_rr, file = paste0(fn_out_path, "df_RR_season_geo_age.tsv"))

# Plot 4 (season) x 3 (geography) grid of heatmaps
df_plot <- df_rr %>%
  mutate(fill_RR = ifelse(N_pairs == 0 | is.na(RR), NA_real_, fill_bound(RR, LB, UB)))

p <- ggplot(df_plot, aes(x = x, y = y, fill = fill_RR)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(season ~ geo_cat) +
  RR_log_grad(LB = LB, UB = UB) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Age-stratified RR by season and geography (Florida vs Mideast)",
       x = "Age group", y = "Age group") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.4, "lines")
  ) +
  coord_equal()

fn_fig_path <- paste0("figs/", scenario, "/season_geo_age/")
dir.create(fn_fig_path, recursive = TRUE, showWarnings = FALSE)
ggsave(paste0(fn_fig_path, "age_RR_season_geo.png"), plot = p, width = 9, height = 11, dpi = 300)
ggsave(paste0(fn_fig_path, "age_RR_season_geo.svg"), plot = p, width = 9, height = 11)

DBI::dbDisconnect(con, shutdown = TRUE)
print("Successfully finished seasonal FL/Mideast age analysis!")
