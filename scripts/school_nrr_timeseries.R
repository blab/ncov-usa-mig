#!/usr/bin/env Rscript
#File: school_nrr_timeseries.R
#Author(s): Amin Bemanian
#Date: 6/9/26
#Description: Normalized RR (nRR) season-by-year line traces for WITHIN-STATE
# identical-sequence pairs pooled across all divisions, using school-relevant age
# bins (0-4, 5-11, 12-17, 18-24, 25-64, 65+). For a chosen index age group, each
# line shows nRR(index, partner_age) over time, normalized to the 25-64/25-64
# within-state diagonal at each time point. A small log-scaled bar panel of total
# within-state pairs sits beneath each trace to convey sample size. The goal is to
# see whether transmission between children changed across the pandemic years.
# One figure per index age group.

library(tidyverse)
library(argparse)
library(duckdb)
library(dbplyr)
library(lubridate)
library(patchwork)

source("scripts/color_schemes.R")
source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")

# Fixed-baseline normalization: nRR_fixed(x,y) = RR(x,y) / RR(ref,ref).
# Applied per stratum so the baseline is taken within the same time point.
normalized_age_rr_fixed <- function(df_rr, baseline_grp){
  rr_baseline <- df_rr %>%
    filter(x == baseline_grp & y == baseline_grp) %>%
    select(RR_baseline = RR, N_baseline = N_pairs) %>%
    mutate(RR_baseline = ifelse(N_baseline == 0, NA, RR_baseline))
  df_rr %>%
    mutate(nRR_fixed = RR / rr_baseline$RR_baseline[1])
}

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--ref_age', type = 'character', default = "25-64", help = 'Reference age bin for fixed-baseline normalization (default 25-64)')
  parser$add_argument('--index_ages', type = 'character', default = "all", help = 'Comma-separated index age bins to make figures for, or "all" (default)')
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE, help = "Whether to exclude possible duplicate pairs, default is FALSE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
ref_age <- args$ref_age
exclude_duplicates <- args$exclude_duplicates

# School-relevant, irregular-width age bins
SCHOOL_LEVELS <- c("0-4", "5-11", "12-17", "18-24", "25-64", "65+")
age_school_bin <- function(a){
  case_when(
    a <= 4  ~ "0-4",
    a <= 11 ~ "5-11",
    a <= 17 ~ "12-17",
    a <= 24 ~ "18-24",
    a <= 64 ~ "25-64",
    TRUE    ~ "65+"
  )
}

if(!(ref_age %in% SCHOOL_LEVELS)){
  stop(sprintf("Reference age '%s' is not one of: %s", ref_age, paste(SCHOOL_LEVELS, collapse = ", ")))
}

assign_season <- function(m){
  case_when(
    m %in% c(3, 4, 5)   ~ "Mar-May",
    m %in% c(6, 7, 8)   ~ "Jun-Aug",
    m %in% c(9, 10, 11) ~ "Sep-Nov",
    TRUE                ~ "Dec-Feb"
  )
}

# Representative date for each season instance. Winter (Dec-Feb) is anchored to
# mid-January, with December rolled into the following year so a single winter
# spans the year boundary cleanly.
season_date <- function(d){
  m <- month(d); y <- year(d)
  case_when(
    m %in% c(3, 4, 5)   ~ make_date(y, 4, 15),
    m %in% c(6, 7, 8)   ~ make_date(y, 7, 15),
    m %in% c(9, 10, 11) ~ make_date(y, 10, 15),
    m == 12             ~ make_date(y + 1, 1, 15),
    TRUE                ~ make_date(y, 1, 15)  # Jan, Feb
  )
}

fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

full_bounds <- c(as.Date("2019-01-01"), as.Date("2025-01-01"))

# Bind division (carries transmit_date through pairs_time); join raw age separately
pairs_div <- con %>%
  bind_pairs_exp("division", time_bounds = full_bounds, exclude_duplicates = exclude_duplicates) %>%
  rename(division.x = x, division.y = y)

age_dict <- tbl(con, "metadata") %>%
  select(strain, age_adj) %>%
  collect()

df_pairs <- pairs_div %>%
  collect() %>%
  left_join(age_dict %>% rename(age.x = age_adj), by = join_by(strain_1 == strain)) %>%
  left_join(age_dict %>% rename(age.y = age_adj), by = join_by(strain_2 == strain)) %>%
  filter(!is.na(age.x) & !is.na(age.y))

DBI::dbDisconnect(con, shutdown = TRUE)

# Restrict to within-state pairs (pooled across all divisions), derive season + bins
df_pairs <- df_pairs %>%
  filter(division.x == division.y) %>%
  mutate(
    season = assign_season(month(transmit_date)),
    season_date = season_date(transmit_date),
    season_label = paste0(season, " ", year(season_date)),
    x = age_school_bin(age.x),
    y = age_school_bin(age.y)
  ) %>%
  filter(!is.na(x) & !is.na(y))

# Compute RR matrix + fixed-baseline nRR per season time point (single geography)
strata <- df_pairs %>% distinct(season_date, season, season_label)

df_nrr <- pmap_dfr(strata, function(season_date, season, season_label){
  subset <- df_pairs %>%
    filter(season_date == !!season_date) %>%
    select(strain_1, strain_2, x, y)
  if(nrow(subset) == 0) return(NULL)
  subset %>%
    calculate_rr_matrix() %>%
    normalized_age_rr_fixed(baseline_grp = ref_age) %>%
    mutate(season_date = !!season_date, season = !!season, season_label = !!season_label)
})

df_nrr <- df_nrr %>%
  mutate(nRR_fixed = ifelse(N_pairs == 0, NA_real_, nRR_fixed)) %>%
  mutate(
    x = factor(x, levels = SCHOOL_LEVELS),
    y = factor(y, levels = SCHOOL_LEVELS)
  )

# Write long table
fn_out_path <- paste0("results/", scenario, "/time_age/")
dir.create(fn_out_path, recursive = TRUE, showWarnings = FALSE)
write_tsv(df_nrr, file = paste0(fn_out_path, "df_nRR_school_season_timeseries.tsv"))

# Total within-state pairs per season (sample size underlying the nRR estimates)
pair_counts <- df_pairs %>%
  group_by(season_date) %>%
  summarise(n_pairs = n(), .groups = "drop")

# Shared x-axis: one labeled break per season time point. Pad limits so edge bars
# (width ~70 days) are not clipped.
x_breaks <- df_nrr %>% distinct(season_date, season_label) %>% arrange(season_date)
date_lims <- range(df_nrr$season_date) + c(-60, 60)

# One figure per index age group
if(args$index_ages == "all"){
  index_ages <- SCHOOL_LEVELS
} else {
  index_ages <- str_split(args$index_ages, ",")[[1]] %>% str_trim()
}

fn_fig_path <- paste0("figs/", scenario, "/age_time/school/")
dir.create(fn_fig_path, recursive = TRUE, showWarnings = FALSE)

for(idx in index_ages){
  df_idx <- df_nrr %>% filter(x == idx, !is.na(nRR_fixed))
  if(nrow(df_idx) == 0){
    message(sprintf("No data for index age %s, skipping.", idx))
    next
  }

  # Pairs that include the index age group (either strain), per season
  idx_counts <- df_pairs %>%
    filter(x == idx | y == idx) %>%
    group_by(season_date) %>%
    summarise(n_index = n(), .groups = "drop")

  p_line <- ggplot(df_idx, aes(x = season_date, y = nRR_fixed, color = y, group = y)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1) +
    scale_color_viridis_d(name = "Partner age", option = "turbo", drop = FALSE) +
    scale_x_date(limits = date_lims, breaks = x_breaks$season_date, labels = x_breaks$season_label) +
    scale_y_log10() +
    labs(
      title = sprintf("Within-state nRR over time for index age %s (ref: %s/%s)", idx, ref_age, ref_age),
      subtitle = "Identical-pair age mixing by season, pooled across all states (within-state pairs)",
      y = "nRR (fixed baseline)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.key.size = unit(0.4, "cm")
    )

  # Total pairs (grey) with the index-involving portion overlaid in red; the red
  # height vs grey top reads as the proportion of pairs including the index age.
  bar_fill_lvls <- c("All within-state pairs", paste0("Includes ", idx))
  p_bar <- ggplot() +
    geom_col(data = pair_counts, aes(x = season_date, y = n_pairs, fill = bar_fill_lvls[1]), width = 70) +
    geom_col(data = idx_counts, aes(x = season_date, y = n_index, fill = bar_fill_lvls[2]), width = 70) +
    scale_fill_manual(name = NULL, values = setNames(c("grey55", "red3"), bar_fill_lvls)) +
    scale_x_date(limits = date_lims, breaks = x_breaks$season_date, labels = x_breaks$season_label) +
    scale_y_log10(breaks = c(1e1, 1e3, 1e5), labels = c("1E1", "1E3", "1E5")) +
    labs(y = "Pairs", x = "Season") +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      legend.key.size = unit(0.4, "cm")
    )

  p <- wrap_plots(list(p_line, p_bar), ncol = 1, heights = c(4, 1), guides = "collect")

  fn_safe <- str_replace_all(idx, "[/+]", "p")
  ggsave(paste0(fn_fig_path, "school_nRR_trace_", fn_safe, ".png"), plot = p, width = 9, height = 7, dpi = 300)
}

cat("School nRR season time series analysis complete.\n")
