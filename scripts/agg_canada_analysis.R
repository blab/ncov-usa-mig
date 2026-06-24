#File: agg_canada_analysis.R
#Author(s): Amin Bemanian
#Date: 6/21/26
#Description: Hypothesis test - calculate geographic RR treating Canada as a SINGLE
#  aggregated geographic unit (instead of splitting it into provinces), alongside
#  all individual US states and Mexico. Goal is to check whether aggregating Canada
#  reveals cross-border RR > 1 signals that were missed (or noisier) when Canada was
#  split province-by-province in the standard state_analysis.R.
#
#  Mechanism: adds an additive, idempotent `agg_geo` column to the metadata table
#  (= "Canada" for any Canadian province, otherwise = division). This column is wiped
#  and rebuilt any time clean_data.R re-writes the metadata table, so it is transient.
#  This lets us reuse the existing, tested bind_pairs_exp / calculate_rr_matrix /
#  calculate_rr_ci pipeline unchanged.
#Arguments:
#--scenario: Scenario corresponding to data files (default CAM_1000)
#--ci: Whether to calculate bootstrap CIs (default TRUE)
#--exclude_duplicates: Whether to exclude possible duplicate pairs (default FALSE)

library(argparse)
library(tidyverse)
library(data.table)
library(duckdb)
library(dbplyr)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', default = "CAM_1000", type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--ci', type = 'logical', default = TRUE, help = "Whether to calculate CIs, default is TRUE")
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE, help = "Whether to exclude possible duplicate pairs, default is FALSE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
ci_flag <- args$ci
exclude_duplicates <- args$exclude_duplicates

fn_db <- paste0("db_files/db_",scenario,".duckdb")
drv <- duckdb(config = list(threads = "15"))
con <- DBI::dbConnect(drv,fn_db)

# --- Build the aggregated-Canada exposure column on the metadata table ---
# Additive + idempotent: "Canada" for any Canadian province, else the division.
DBI::dbExecute(con, "ALTER TABLE metadata ADD COLUMN IF NOT EXISTS agg_geo VARCHAR")
DBI::dbExecute(con, "UPDATE metadata SET agg_geo = CASE WHEN country = 'Canada' THEN 'Canada' ELSE division END")
print("Added/updated agg_geo column (Canada aggregated into a single unit).")

# --- RR matrix using aggregated geography ---
agg_rr <- con %>%
  bind_pairs_exp("agg_geo", exclude_duplicates = exclude_duplicates) %>%
  calculate_rr_matrix() %>%
  collect()

if(ci_flag){
  agg_rr_ci <- calculate_rr_ci(con, "agg_geo", exclude_duplicates = exclude_duplicates)
  agg_rr <- inner_join(agg_rr, agg_rr_ci, by=join_by(x,y))
}

fn_rr <- paste0("results/",scenario,"/df_RR_agg_canada.tsv")
print(fn_rr)
readr::write_tsv(agg_rr, file=fn_rr)

# --- Hypothesis summary: Canada <-> US cross-border RR ---
# Cross-border = one deme is aggregated "Canada", the other a US state (i.e. not
# Canada and not Mexico). Each unordered pair appears twice (x,y) and (y,x); we keep
# the rows where x == "Canada" to get one row per US state.
canada_cross <- agg_rr %>%
  filter(x == "Canada", y != "Canada", y != "Mexico") %>%
  arrange(desc(RR))

ci_sig <- "ci_lb" %in% names(canada_cross)
n_states <- nrow(canada_cross)
n_rr_gt1 <- canada_cross %>% filter(RR > 1) %>% nrow()

cat("\n===== AGGREGATED-CANADA CROSS-BORDER RR SUMMARY =====\n")
cat(sprintf("US states/territories compared against aggregated Canada: %d\n", n_states))
cat(sprintf("  with point-estimate RR > 1: %d\n", n_rr_gt1))
if(ci_sig){
  n_sig_gt1 <- canada_cross %>% filter(ci_lb > 1) %>% nrow()
  cat(sprintf("  with RR > 1 AND CI lower bound > 1 (significant): %d\n", n_sig_gt1))
}

cat("\nTop cross-border RR values (aggregated Canada):\n")
if(ci_sig){
  print(canada_cross %>% select(x, y, RR, ci_lb, ci_ub) %>% head(15))
}else{
  print(canada_cross %>% select(x, y, RR) %>% head(15))
}

# --- Contrast against province-level results (if available) ---
# Loads the standard state_analysis.R output and reports how the Canada<->US signal
# looked when Canada was split into provinces, for direct comparison.
fn_prov <- paste0("results/",scenario,"/df_RR_by_state.tsv")
if(file.exists(fn_prov)){
  canada_provinces <- con %>%
    tbl("metadata") %>%
    filter(country == "Canada") %>%
    distinct(division) %>%
    pull(division)

  prov_rr <- readr::read_tsv(fn_prov, show_col_types = FALSE)
  prov_cross <- prov_rr %>%
    filter(x %in% canada_provinces, !(y %in% canada_provinces), y != "Mexico")

  prov_ci <- "ci_lb" %in% names(prov_cross)
  cat("\n===== PROVINCE-LEVEL CROSS-BORDER RR (for contrast) =====\n")
  cat(sprintf("Canadian province <-> US state pairs: %d\n", nrow(prov_cross)))
  cat(sprintf("  with point-estimate RR > 1: %d\n", prov_cross %>% filter(RR > 1) %>% nrow()))
  if(prov_ci){
    cat(sprintf("  with RR > 1 AND CI lower bound > 1 (significant): %d\n",
                prov_cross %>% filter(ci_lb > 1) %>% nrow()))
  }
  cat(sprintf("  mean province<->US RR: %.3f | mean aggregated-Canada<->US RR: %.3f\n",
              mean(prov_cross$RR, na.rm=TRUE), mean(canada_cross$RR, na.rm=TRUE)))
}else{
  cat(sprintf("\n(Province-level file %s not found - skipping contrast.)\n", fn_prov))
}

DBI::dbDisconnect(con, shutdown=TRUE)
print("Successfully finished aggregated Canada analysis!")
