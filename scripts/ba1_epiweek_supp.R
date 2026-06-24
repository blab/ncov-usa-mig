#File: ba1_epiweek_supp.R
#Author(s): Amin Bemanian
#Date: 6/21/26
#Description: Supplementary figure documenting the BA.1 (Omicron 21K) sampling-bias
#  artifact behind the Oklahoma<->Mexico identical-pair RR spike. Plots the temporal
#  distribution of sequencing effort (percent of each geography's total sequences) by
#  epidemiological week for Oklahoma, Texas, California, and Mexico. Weeks falling in
#  the BA.1 wave are highlighted in red; facet titles report each geography's total N.
#  Shows that Oklahoma packed a disproportionate share of its sequencing into the BA.1
#  window relative to the genuine border states.
#Arguments:
#--scenario: Scenario corresponding to data files (default CAM_1000)

library(argparse)
library(tidyverse)
library(duckdb)
library(dbplyr)
library(scales)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', default = "CAM_1000", type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# BA.1 (Omicron 21K) wave window - the low-diversity global sweep of winter 2021-22
BA1_START <- as.Date("2021-12-01")
BA1_END   <- as.Date("2022-03-31")

# Geographies of interest, in display order
GEOS <- c("Oklahoma", "Texas", "California", "Mexico")

BA1_RED  <- "#cb181d"
OTHER_GREY <- "grey70"

fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

meta <- tbl(con, "metadata") %>%
  filter(division %in% GEOS) %>%
  select(strain, division, date) %>%
  collect()

DBI::dbDisconnect(con, shutdown = TRUE)

# Parse dates; drop partial-precision records (year-only / year-month) that cannot be
# assigned to an epi week. These are <1% of sequences for each geography.
seq_wk <- meta %>%
  mutate(date = ymd(date, quiet = TRUE)) %>%
  filter(!is.na(date)) %>%
  # Epidemiological week = MMWR week, which starts on Sunday
  mutate(epiweek = floor_date(date, unit = "week", week_start = 7)) %>%
  count(division, epiweek, name = "n_seq") %>%
  group_by(division) %>%
  mutate(pct = 100 * n_seq / sum(n_seq)) %>%
  ungroup() %>%
  mutate(ba1 = epiweek >= BA1_START & epiweek <= BA1_END,
         ba1 = factor(if_else(ba1, "BA.1 window (21K)", "Other"),
                      levels = c("Other", "BA.1 window (21K)")))

# Facet labels with total (plotted) N per geography, ordered as requested
facet_n <- seq_wk %>%
  group_by(division) %>%
  summarize(total_n = sum(n_seq), .groups = "drop") %>%
  mutate(facet_label = paste0(division, " (n = ", comma(total_n), ")"))

seq_wk <- seq_wk %>%
  left_join(facet_n, by = "division") %>%
  mutate(facet_label = factor(facet_label,
                              levels = facet_n$facet_label[match(GEOS, facet_n$division)]))

p_ba1 <- ggplot(seq_wk, aes(x = epiweek, y = pct, fill = ba1)) +
  # Faint band marking the BA.1 window across all facets
  annotate("rect", xmin = BA1_START, xmax = BA1_END, ymin = -Inf, ymax = Inf,
           fill = BA1_RED, alpha = 0.08) +
  geom_col(width = 7) +
  facet_wrap(~ facet_label, ncol = 1, scales = "fixed") +
  scale_fill_manual(values = c("Other" = OTHER_GREY, "BA.1 window (21K)" = BA1_RED),
                    name = NULL) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y", name = "Epidemiological week") +
  scale_y_continuous(name = "Percent of geography's total sequences (%)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 8),
        strip.text = element_text(face = "bold"))

FN_PATH <- paste0("figs/", scenario, "/supp/")
dir.create(FN_PATH, recursive = TRUE, showWarnings = FALSE)

ggsave(paste0(FN_PATH, "ba1_epiweek_seqdist.jpg"),
       p_ba1, height = 9, width = 7, units = "in", dpi = 300)

fn_supp_pdf <- "manuscript/figures/supp/ba1_epiweek_seqdist.pdf"
ggsave(fn_supp_pdf, p_ba1, height = 9, width = 7, units = "in")

print(paste0("Saved: ", FN_PATH, "ba1_epiweek_seqdist.jpg and ", fn_supp_pdf))
print("Successfully finished BA.1 epiweek supplement figure!")
