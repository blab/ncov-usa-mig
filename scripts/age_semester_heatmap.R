# Age RR / nRR heatmaps with 3-year age bins, faceted by school semester
# from Fall 2020 through Summer 2024 (early 2020 months dropped due to low
# sequence counts). Semesters: Spring (Jan 1 - May 31), Summer break
# (Jun 1 - Jul 31), and Fall (Aug 1 - Dec 31) of each year.
# Produces two figures:
#   1) raw RR
#   2) nRR_fixed normalized to the 39-41y within-group RR (per semester)

library(argparse)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
library(duckdb)
library(dbplyr)

source("scripts/calculate_rr_matrix.R")
source("scripts/color_schemes.R")

# nRR with a fixed baseline group: nRR_fixed = RR(x,y) / RR(baseline,baseline),
# normalized within each grouping level (here: per semester).
normalized_age_rr_fixed <- function(df_rr, baseline_grp, group_vars){
  rr_baseline <- df_rr %>%
    filter(x == baseline_grp & y == baseline_grp) %>%
    select(all_of(group_vars), RR_baseline = RR, N_baseline = N_pairs) %>%
    mutate(RR_baseline = ifelse(N_baseline == 0, NA, RR_baseline))
  df_rr %>%
    left_join(rr_baseline, by = group_vars) %>%
    mutate(nRR_fixed = RR / RR_baseline) %>%
    select(-RR_baseline, -N_baseline)
}

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

BIN_WIDTH <- 3           # Age bin width in years (0-2, 3-5, ...)
AGE_CAP <- 80            # Cap at age 80
BASELINE <- 39           # nRR baseline bin lower bound -> 39-41y
LB <- 0.5; UB <- 2.0     # RR/nRR display bounds (broad to handle per-semester noise)
MIN_PAIRS <- 5           # Block out (gray) cells with fewer identical pairs than this

# 3-year-bin label from the bin's lower bound
bin_label <- function(lo) sprintf("%02d-%02dy", lo, lo + BIN_WIDTH - 1)

# School semesters from Fall 2020 through Summer 2024 (chronological).
# Built as season x year, then trimmed to the Fall 2020 - Summer 2024 window.
season_defs <- tibble::tribble(
  ~season,   ~mstart,   ~mend,
  "Spring",  "01-01",   "05-31",
  "Summer",  "06-01",   "07-31",
  "Fall",    "08-01",   "12-31"
)
# Columns ordered by academic-year progression (Fall -> Spring -> Summer) so the
# Fall 2020 - Summer 2024 window fills a complete grid with no blank panels.
SEASON_LEVELS <- c("Fall", "Spring", "Summer")
semesters <- tidyr::expand_grid(year = 2020:2024, season_defs) %>%
  mutate(
    season = factor(season, levels = SEASON_LEVELS),
    start  = paste0(year, "-", mstart),
    end    = paste0(year, "-", mend),
    label  = paste(season, year),
    # Academic year: Fall belongs to its own calendar year; the following
    # Spring/Summer belong to the academic year that started the prior Fall.
    acad_start = ifelse(season == "Fall", year, year - 1L),
    acad_year  = paste0(acad_start, "-", substr(as.character(acad_start + 1L), 3, 4))
  ) %>%
  filter(!(year == 2020 & season %in% c("Spring", "Summer"))) %>%  # drop low-count early 2020
  filter(!(year == 2024 & season == "Fall")) %>%                   # end at Summer 2024
  arrange(acad_start, season)

fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(config = list(threads = "15")), fn_db)

# Bind pairs to 3-year age bins (keyed by bin lower bound) for a time window.
# Mirrors bind_pairs_exp but injects the binning expression on age_adj.
bind_pairs_age_bin <- function(con, time_bounds){
  time_lb <- time_bounds[1]; time_ub <- time_bounds[2]
  exp_dict <- tbl(con, "metadata") %>%
    filter(!is.na(age_adj), age_adj >= 0, age_adj <= AGE_CAP) %>%
    transmute(strain, exp = (age_adj %/% BIN_WIDTH) * BIN_WIDTH)
  tbl(con, "pairs_time") %>%
    filter(transmit_date > time_lb, transmit_date < time_ub) %>%
    inner_join(exp_dict, join_by(strain_1 == strain)) %>% rename(x = exp) %>%
    inner_join(exp_dict, join_by(strain_2 == strain)) %>% rename(y = exp) %>%
    filter(!is.na(x), !is.na(y))
}

# Compute the 3-year-bin age RR for each semester window
df_all <- map_dfr(seq_len(nrow(semesters)), function(i){
  message("Calculating ", semesters$label[i])
  tb <- c(as.Date(semesters$start[i]) - 1, as.Date(semesters$end[i]) + 1)
  bind_pairs_age_bin(con, tb) %>%
    calculate_rr_matrix() %>%
    collect() %>%
    mutate(semester = semesters$label[i])
})

# Weekly sequence counts (valid dates only) for the per-academic-year line plots
weekly <- DBI::dbGetQuery(con, "
  SELECT date_trunc('week', TRY_CAST(date AS DATE)) AS wk, COUNT(*) AS n
  FROM metadata
  WHERE TRY_CAST(date AS DATE) IS NOT NULL
  GROUP BY 1 ORDER BY 1")

DBI::dbDisconnect(con, shutdown = TRUE)

# Assign each week to a season / academic year (matching the heatmap windows)
weekly <- weekly %>%
  mutate(
    yr  = as.integer(format(wk, "%Y")),
    mon = as.integer(format(wk, "%m")),
    season = case_when(mon >= 8 ~ "Fall", mon <= 5 ~ "Spring", TRUE ~ "Summer"),
    acad_start = ifelse(season == "Fall", yr, yr - 1L),
    acad_year  = paste0(acad_start, "-", substr(as.character(acad_start + 1L), 3, 4))
  ) %>%
  filter(acad_start %in% 2020:2023) %>%
  mutate(season = factor(season, levels = SEASON_LEVELS))

# Add baseline-normalized nRR (per semester)
df_all <- normalized_age_rr_fixed(df_all, baseline_grp = BASELINE, group_vars = "semester")

# Order bins / facets and mask under-sampled cells
bin_levels <- sort(unique(c(df_all$x, df_all$y)))
df_all <- df_all %>%
  left_join(semesters %>% select(label, acad_year, season), by = c("semester" = "label")) %>%
  mutate(
    x_lab = factor(bin_label(x), levels = bin_label(bin_levels)),
    y_lab = factor(bin_label(y), levels = bin_label(bin_levels)),
    season = factor(season, levels = SEASON_LEVELS),
    # Block out under-sampled cells: NA fill renders as gray (na.value)
    fill_RR  = ifelse(N_pairs >= MIN_PAIRS, fill_bound(RR, LB, UB), NA_real_),
    fill_nRR = ifelse(N_pairs >= MIN_PAIRS, fill_bound(nRR_fixed, LB, UB), NA_real_)
  )

# Axis breaks every 15 years to avoid clutter
age_breaks <- bin_label(bin_levels[bin_levels %% 15 == 0])

make_heatmap <- function(fill_col, title, legend_name){
  ggplot(df_all, aes(x = x_lab, y = y_lab, fill = .data[[fill_col]])) +
    geom_tile() +
    facet_grid(acad_year ~ season) +
    RR_log_grad(LB = LB, UB = UB) +
    scale_x_discrete(breaks = age_breaks) +
    scale_y_discrete(breaks = age_breaks) +
    coord_equal() +
    labs(title = title, x = "Age", y = "Age", fill = legend_name) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    )
}

out_dir <- paste0("figs/", scenario, "/age_heatmaps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fn_rr <- paste0(out_dir, "/age_semester_heatmap.png")
ggsave(fn_rr, make_heatmap("fill_RR", "Age RR by school semester (3-year bins, Fall 2020 - Summer 2024)", "RR"),
       width = 9, height = 12, dpi = 300)
message("Saved ", fn_rr)

fn_nrr <- paste0(out_dir, "/age_semester_heatmap_nRR.png")
ggsave(fn_nrr, make_heatmap("fill_nRR", "Age nRR (baseline 39-41y) by school semester (3-year bins, Fall 2020 - Summer 2024)", "nRR"),
       width = 9, height = 12, dpi = 300)
message("Saved ", fn_nrr)

# --- Composite: per-academic-year nRR heatmap row + weekly sequence-count line ---
ay_levels <- unique(semesters$acad_year)  # chronological

# One academic-year row: 3 season heatmaps over a 3-season weekly-count line plot
hm_row <- function(ay, fill_col){
  ggplot(filter(df_all, acad_year == ay), aes(x_lab, y_lab, fill = .data[[fill_col]])) +
    geom_tile() +
    facet_wrap(~ season, nrow = 1, drop = FALSE) +
    RR_log_grad(LB = LB, UB = UB) +
    scale_x_discrete(breaks = age_breaks) +
    scale_y_discrete(breaks = age_breaks) +
    coord_equal() +
    labs(title = ay, x = NULL, y = "Age") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "none"
    )
}

ln_row <- function(ay){
  ggplot(filter(weekly, acad_year == ay), aes(wk, n)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ season, nrow = 1, scales = "free_x", drop = FALSE) +
    # Shared log10 y-axis across every panel/block (ceiling 1E6 to keep the
    # ~129k peak weeks; floor 1E3, below the in-window weekly minimum of ~2k)
    scale_y_log10(
      limits = c(1e3, 1e6),
      breaks = 10^(3:6),
      labels = c("1E3", "1E4", "1E5", "1E6")
    ) +
    labs(x = NULL, y = "Seq/wk") +
    theme_bw() +
    theme(
      strip.text = element_blank(),
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 4.5)
    )
}

block <- function(ay, fill_col){
  (hm_row(ay, fill_col) / ln_row(ay)) + plot_layout(heights = c(5, 1.4))
}

composite <- wrap_plots(lapply(ay_levels, block, fill_col = "fill_nRR"), ncol = 1) +
  plot_annotation(
    title = "Age nRR (baseline 39-41y) by school semester with weekly sequence counts",
    theme = theme(plot.title = element_text(face = "bold"))
  )

fn_comp <- paste0(out_dir, "/age_semester_heatmap_nRR_counts.png")
ggsave(fn_comp, composite, width = 9, height = 16, dpi = 200)
message("Saved ", fn_comp)
