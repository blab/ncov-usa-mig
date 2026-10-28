#File: age_season_nRR_traces.R
#Author(s): Amin Bemanian
#Date: 6/9/26
#Description: Normalized RR (nRR) time-series line plots for identical-sequence
# pairs restricted to Florida and the Mideast (BEA region), stratified by
# geographic pairing (FL/FL, Mideast/Mideast, FL/Mideast). Unlike the pooled
# heatmap analysis, time is resolved by season-within-year so we can trace how
# age-mixing evolves over the pandemic. For a chosen index age group, each line
# shows nRR(index, partner_age) over time for a partner age group. nRR is the
# fixed-baseline normalization against the 30-34y/30-34y diagonal within the
# same geography and same time point. One figure per index age group, faceted
# into 3 rows by geography.

library(argparse)
library(tidyverse)
library(duckdb)
library(dbplyr)
library(lubridate)
library(patchwork)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/color_schemes.R")

# Fixed-baseline normalization (matches normalized_age_rr_fixed in
# age_time_RR_analysis.R; inlined to avoid that script's top-level side effects).
# Applied per stratum so the baseline RR(ref,ref) is taken within the same
# geography and time point.
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
  parser$add_argument('--bin_size', type = 'integer', default = 5, help = 'Age bin width in years (default 5)')
  parser$add_argument('--age_cap', type = 'integer', default = 85, help = 'Top-coded age; ages >= cap pooled into a single bin (default 85)')
  parser$add_argument('--ref_age', type = 'character', default = "30-34y", help = 'Reference age bin for fixed-baseline normalization (default 30-34y)')
  parser$add_argument('--index_ages', type = 'character', default = "all", help = 'Comma-separated index age bins to make figures for, or "all" (default)')
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE, help = "Whether to exclude possible duplicate pairs, default is FALSE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
bin_size <- args$bin_size
age_cap <- args$age_cap
ref_age <- args$ref_age
exclude_duplicates <- args$exclude_duplicates

# Configurable, top-coded age binning
age_bin <- function(a, bin = 5, cap = 85){
  ifelse(a >= cap,
         sprintf("%d+y", cap),
         sprintf("%02d-%02dy", (a %/% bin) * bin, (a %/% bin) * bin + (bin - 1)))
}
age_bin_levels <- function(bin = 5, cap = 85){
  starts <- seq(0, cap - bin, by = bin)
  c(sprintf("%02d-%02dy", starts, starts + (bin - 1)), sprintf("%d+y", cap))
}

GEO_LEVELS <- c("FL/FL", "Mideast/Mideast", "FL/Mideast")

assign_season <- function(m){
  case_when(
    m %in% c(3, 4, 5)   ~ "Mar-May",
    m %in% c(6, 7, 8)   ~ "Jun-Aug",
    m %in% c(9, 10, 11) ~ "Sep-Nov",
    TRUE                ~ "Dec-Feb"
  )
}

# Representative date for each season instance. Winter (Dec-Feb) is anchored to
# mid-January, with December rolled forward into the following year so a single
# winter spans the year boundary cleanly.
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

# Geography binds (pairs_div carries transmit_date through pairs_time)
pairs_div <- con %>%
  bind_pairs_exp("division", time_bounds = full_bounds, exclude_duplicates = exclude_duplicates) %>%
  rename(division.x = x, division.y = y)
pairs_reg <- con %>%
  bind_pairs_exp("bea_reg", time_bounds = full_bounds, exclude_duplicates = exclude_duplicates) %>%
  rename(bea_reg.x = x, bea_reg.y = y)

age_dict <- tbl(con, "metadata") %>%
  select(strain, age_adj) %>%
  collect()

df_pairs <- pairs_div %>%
  left_join(pairs_reg %>% select(strain_1, strain_2, bea_reg.x, bea_reg.y), by = join_by(strain_1, strain_2)) %>%
  collect() %>%
  left_join(age_dict %>% rename(age.x = age_adj), by = join_by(strain_1 == strain)) %>%
  left_join(age_dict %>% rename(age.y = age_adj), by = join_by(strain_2 == strain)) %>%
  filter(!is.na(age.x) & !is.na(age.y))

DBI::dbDisconnect(con, shutdown = TRUE)

# Classify geography (restrict to FL + Mideast strains)
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
    season = assign_season(month(transmit_date)),
    season_date = season_date(transmit_date)
  )

# Age binning
bin_levels <- age_bin_levels(bin_size, age_cap)
df_pairs <- df_pairs %>%
  mutate(
    x = age_bin(age.x, bin_size, age_cap),
    y = age_bin(age.y, bin_size, age_cap)
  ) %>%
  filter(!is.na(x) & !is.na(y))

if(!(ref_age %in% bin_levels)){
  stop(sprintf("Reference age '%s' is not a valid bin for bin_size=%d, age_cap=%d. Valid bins: %s",
               ref_age, bin_size, age_cap, paste(bin_levels, collapse = ", ")))
}

# Compute RR matrix + fixed-baseline nRR per (season_date x geo_cat) stratum
strata <- df_pairs %>% distinct(season_date, season, geo_cat)

df_nrr <- pmap_dfr(strata, function(season_date, season, geo_cat){
  subset <- df_pairs %>%
    filter(season_date == !!season_date & geo_cat == !!geo_cat) %>%
    select(strain_1, strain_2, x, y)
  if(nrow(subset) == 0) return(NULL)
  subset %>%
    calculate_rr_matrix() %>%
    normalized_age_rr_fixed(baseline_grp = ref_age) %>%
    mutate(season_date = !!season_date, season = !!season, geo_cat = !!geo_cat)
})

# Drop cells with no observed pairs (RR uses +1 smoothing; those nRR are noise)
df_nrr <- df_nrr %>%
  mutate(nRR_fixed = ifelse(N_pairs == 0, NA_real_, nRR_fixed)) %>%
  mutate(
    geo_cat = factor(geo_cat, levels = GEO_LEVELS),
    x = factor(x, levels = bin_levels),
    y = factor(y, levels = bin_levels)
  )

# Write long table
fn_out_path <- paste0("results/", scenario, "/season_geo_age/")
dir.create(fn_out_path, recursive = TRUE, showWarnings = FALSE)
write_tsv(df_nrr, file = paste0(fn_out_path, "df_nRR_traces_season_geo_age.tsv"))

# Total identical pairs per geography x season (sample size underlying the nRR
# estimates in each panel; index-age independent)
pair_counts <- df_pairs %>%
  group_by(geo_cat, season_date) %>%
  summarise(n_pairs = n(), .groups = "drop")

# Pad limits so edge bars (width ~70 days) are not clipped
date_lims <- range(df_nrr$season_date) + c(-60, 60)

# One figure per index age group: trace nRR(index, partner) over time, 3 rows by
# geo, each with a small bar panel of total pairs underneath to convey sample size
if(args$index_ages == "all"){
  index_ages <- bin_levels
} else {
  index_ages <- str_split(args$index_ages, ",")[[1]] %>% str_trim()
}

fn_fig_path <- paste0("figs/", scenario, "/season_geo_age/nRR_traces/")
dir.create(fn_fig_path, recursive = TRUE, showWarnings = FALSE)

# Build the line panel + bar panel pair for one geography
make_geo_block <- function(df_line, df_bar, geo, show_x){
  p_line <- ggplot(df_line, aes(x = season_date, y = nRR_fixed, color = y, group = y)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.8) +
    scale_color_viridis_d(name = "Partner age", option = "turbo", drop = FALSE) +
    scale_x_date(limits = date_lims, date_breaks = "1 year", date_labels = "%Y") +
    scale_y_log10() +
    labs(title = geo, y = "nRR") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  p_bar <- ggplot(df_bar, aes(x = season_date, y = n_pairs)) +
    geom_col(fill = "grey45", width = 70) +
    scale_x_date(limits = date_lims, date_breaks = "1 year", date_labels = "%Y") +
    scale_y_log10(
      breaks = c(1e1, 1e3, 1e5),
      labels = c("1E1", "1E3", "1E5")
    ) +
    labs(y = "Pairs", x = if(show_x) "Season" else NULL) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 8),
      axis.text.x = if(show_x) element_text() else element_blank(),
      axis.ticks.x = if(show_x) element_line() else element_blank()
    )

  list(line = p_line, bar = p_bar)
}

for(idx in index_ages){
  df_idx <- df_nrr %>% filter(x == idx, !is.na(nRR_fixed))
  if(nrow(df_idx) == 0){
    message(sprintf("No data for index age %s, skipping.", idx))
    next
  }

  panels <- list()
  for(i in seq_along(GEO_LEVELS)){
    geo <- GEO_LEVELS[i]
    blk <- make_geo_block(
      df_line = df_idx %>% filter(geo_cat == geo),
      df_bar  = pair_counts %>% filter(geo_cat == geo),
      geo = geo,
      show_x = (i == length(GEO_LEVELS))  # only label x on the bottom-most bar
    )
    panels <- c(panels, list(blk$line, blk$bar))
  }

  p <- wrap_plots(panels, ncol = 1, heights = rep(c(4, 1), length(GEO_LEVELS)), guides = "collect") +
    plot_annotation(
      title = sprintf("nRR over time for index age %s (ref: %s/%s within geo)", idx, ref_age, ref_age),
      subtitle = "Identical-pair age mixing by season; bars show total pairs per geography (Florida vs Mideast)"
    ) &
    theme(legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 6))

  fn_safe <- str_replace_all(idx, "[/+]", "p")
  ggsave(paste0(fn_fig_path, "age_nRR_trace_", fn_safe, ".png"), plot = p, width = 8, height = 10, dpi = 300)
}

print("Successfully finished seasonal FL/Mideast nRR trace analysis!")
