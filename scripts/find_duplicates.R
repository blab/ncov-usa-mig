#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(duckdb)
library(dbplyr)
library(readr)
library(lubridate)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = 'CAM_1000', help = 'Scenario name')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

# Query for potential duplicate sequences (same host)
# Criteria: Pair exists in pairs table with distance=0, AND metadata matches:
#   - same age_adj, same division, within 28 days
#   - sex must match OR both be NA
#   - location must match OR both be NA
df_duplicates <- tbl(con, "pairs") %>%
  # Join metadata for strain_1
  left_join(
    tbl(con, "metadata") %>%
      filter(str_detect(date, "^\\d{4}-\\d{2}-\\d{2}$")) %>%
      mutate(date_1 = sql("CAST(date AS DATE)")) %>%
      select(strain_1 = strain, date_1, division_1 = division, location_1 = location,
             age_adj_1 = age_adj, sex_1 = sex),
    by = "strain_1"
  ) %>%
  # Join metadata for strain_2
  left_join(
    tbl(con, "metadata") %>%
      filter(str_detect(date, "^\\d{4}-\\d{2}-\\d{2}$")) %>%
      mutate(date_2 = sql("CAST(date AS DATE)")) %>%
      select(strain_2 = strain, date_2, division_2 = division, location_2 = location,
             age_adj_2 = age_adj, sex_2 = sex),
    by = "strain_2"
  ) %>%
  # Apply filters
  filter(
    n_mutations == 0,                           # Identical sequences
    division_1 == division_2,                   # Same division (state/province)
    age_adj_1 == age_adj_2,                     # Same age
    abs(date_1 - date_2) <= 28,                 # Within 28 days
    (sex_1 == sex_2) | (is.na(sex_1) & is.na(sex_2)),           # Sex matches or both NA
    (location_1 == location_2) | (is.na(location_1) & is.na(location_2))  # Location (sub-state/province) matches or both NA
  ) %>%
  collect()

cat("Total potential duplicate pairs found:", nrow(df_duplicates), "\n")

# Categorize by metadata availability
df_duplicates <- df_duplicates %>%
  mutate(
    both_sex_valid = !is.na(sex_1) & !is.na(sex_2),
    both_location_valid = !is.na(location_1) & !is.na(location_2),
    year = lubridate::year(pmin(date_1, date_2, na.rm = TRUE))
  )

cat("\nPairs with both sex values valid:", sum(df_duplicates$both_sex_valid), "\n")
cat("Pairs with both location values valid:", sum(df_duplicates$both_location_valid), "\n")
cat("Pairs with both sex AND location valid:", sum(df_duplicates$both_sex_valid & df_duplicates$both_location_valid), "\n")

cat("\nBreakdown by year (of earliest sequence):\n")
df_duplicates %>%
  group_by(year) %>%
  summarise(n_duplicate_pairs = n()) %>%
  arrange(year) %>%
  print()

# Summary statistics
cat("\nSummary by division:\n")
df_duplicates %>%
  group_by(division_1) %>%
  summarise(n_duplicate_pairs = n()) %>%
  arrange(desc(n_duplicate_pairs)) %>%
  print(n = 20)

cat("\nSummary by location (top 20):\n")
df_duplicates %>%
  filter(both_location_valid) %>%
  group_by(location_1) %>%
  summarise(n_duplicate_pairs = n()) %>%
  arrange(desc(n_duplicate_pairs)) %>%
  print(n = 20)

cat("\nSummary by age:\n")
df_duplicates %>%
  group_by(age_adj_1) %>%
  summarise(n_duplicate_pairs = n()) %>%
  arrange(desc(n_duplicate_pairs)) %>%
  print(n = 20)

cat("\nSummary by sex (when both valid):\n")
df_duplicates %>%
  filter(both_sex_valid) %>%
  group_by(sex_1) %>%
  summarise(n_duplicate_pairs = n()) %>%
  arrange(desc(n_duplicate_pairs)) %>%
  print()

cat("\nSummary by date difference:\n")
df_duplicates %>%
  mutate(date_diff = abs(date_1 - date_2)) %>%
  group_by(date_diff) %>%
  summarise(n_duplicate_pairs = n()) %>%
  arrange(date_diff) %>%
  print(n = 29)

# Save results
out_dir <- paste0("results/", scenario, "/duplicates/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(df_duplicates, paste0(out_dir, "potential_duplicate_pairs.tsv"))
cat("\nResults saved to:", paste0(out_dir, "potential_duplicate_pairs.tsv"), "\n")

# Update pairs table with possible_duplicates column
cat("\nUpdating pairs table with possible_duplicates column...\n")

# First, add the column and set all to FALSE
DBI::dbExecute(con, "ALTER TABLE pairs ADD COLUMN IF NOT EXISTS possible_duplicates BOOLEAN DEFAULT FALSE")
DBI::dbExecute(con, "UPDATE pairs SET possible_duplicates = FALSE")

# Create a temporary table with the duplicate pairs
duplicate_pairs <- df_duplicates %>%
  select(strain_1, strain_2) %>%
  distinct()

dbWriteTable(con, "temp_duplicates", duplicate_pairs, temporary = TRUE, overwrite = TRUE)

# Update pairs table to mark duplicates as TRUE
n_updated <- DBI::dbExecute(con,
  "UPDATE pairs
   SET possible_duplicates = TRUE
   WHERE EXISTS (
     SELECT 1 FROM temp_duplicates
     WHERE temp_duplicates.strain_1 = pairs.strain_1
       AND temp_duplicates.strain_2 = pairs.strain_2
   )")

cat("Updated", n_updated, "pairs as possible duplicates in pairs table\n")

# If pairs_time table exists, update it as well
if("pairs_time" %in% DBI::dbListTables(con)){
  cat("\nUpdating pairs_time table with possible_duplicates column...\n")

  # Add column if it doesn't exist and set all to FALSE
  DBI::dbExecute(con, "ALTER TABLE pairs_time ADD COLUMN IF NOT EXISTS possible_duplicates BOOLEAN DEFAULT FALSE")
  DBI::dbExecute(con, "UPDATE pairs_time SET possible_duplicates = FALSE")

  # Update pairs_time table to mark duplicates as TRUE
  n_updated_time <- DBI::dbExecute(con,
    "UPDATE pairs_time
     SET possible_duplicates = TRUE
     WHERE EXISTS (
       SELECT 1 FROM temp_duplicates
       WHERE temp_duplicates.strain_1 = pairs_time.strain_1
         AND temp_duplicates.strain_2 = pairs_time.strain_2
     )")

  cat("Updated", n_updated_time, "pairs as possible duplicates in pairs_time table\n")
} else {
  cat("\nNote: pairs_time table does not exist yet. Run time_prep.R to create it.\n")
  cat("The possible_duplicates column will be included when pairs_time is created from pairs table.\n")
}

DBI::dbDisconnect(con, shutdown=TRUE)
