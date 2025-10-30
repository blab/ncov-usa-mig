library(tidyverse)
library(dbplyr)
library(duckdb)
library(tictoc)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")

# Calculate normalized RR (nRR)
# nRR = RR(x,y) / mean(RR(x,x), RR(y,y))
normalized_age_rr_diag <- function(df_rr, group_vars = NULL){
  # First: filter to the diagonal (RR(x,x)) entries
  if(is.null(group_vars)){
    rr_diag <- df_rr %>%
      filter(x == y) %>%
      select(age = x, RR_diag = RR, N_diag = N_pairs) %>%
      mutate(RR_diag = ifelse(N_diag == 0, NA, RR_diag))

    # Join to get RR(x,x) and RR(y,y) for each row
    df_rr <- df_rr %>%
      left_join(rr_diag, by = c("x" = "age")) %>%
      rename(RR_xx = RR_diag) %>%
      left_join(rr_diag, by = c("y" = "age")) %>%
      rename(RR_yy = RR_diag) %>%
      mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
  } else {
    rr_diag <- df_rr %>%
      filter(x == y) %>%
      select(all_of(group_vars), age = x, RR_diag = RR, N_diag = N_pairs) %>%
      mutate(RR_diag = ifelse(N_diag == 0, NA, RR_diag))

    # Join to get RR(x,x) and RR(y,y) for each row
    df_rr <- df_rr %>%
      left_join(rr_diag, by = c(group_vars, "x" = "age")) %>%
      rename(RR_xx = RR_diag) %>%
      left_join(rr_diag, by = c(group_vars, "y" = "age")) %>%
      rename(RR_yy = RR_diag) %>%
      mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
  }

  return(df_rr)
}

# Calculate normalized RR (nRR) with fixed baseline
# nRR_fixed = RR(x,y) / RR(baseline, baseline)
normalized_age_rr_fixed <- function(df_rr, baseline_grp, group_vars = NULL){
  # Get the baseline diagonal RR(baseline, baseline)
  if(is.null(group_vars)){
    rr_baseline <- df_rr %>%
      filter(x == baseline_grp & y == baseline_grp) %>%
      select(RR_baseline = RR, N_baseline = N_pairs) %>%
      mutate(RR_baseline = ifelse(N_baseline == 0, NA, RR_baseline))

    # Normalize all RR values by the baseline
    df_rr <- df_rr %>%
      mutate(nRR_fixed = RR / rr_baseline$RR_baseline[1])
  } else {
    rr_baseline <- df_rr %>%
      filter(x == baseline_grp & y == baseline_grp) %>%
      select(all_of(group_vars), RR_baseline = RR, N_baseline = N_pairs) %>%
      mutate(RR_baseline = ifelse(N_baseline == 0, NA, RR_baseline))

    # Join baseline RR and normalize
    df_rr <- df_rr %>%
      left_join(rr_baseline, by = group_vars) %>%
      mutate(nRR_fixed = RR / RR_baseline) %>%
      select(-RR_baseline, -N_baseline)
  }

  return(df_rr)
}

collect_args <- function(){
  parser <- argparse::ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--ci', type = 'logical', default = TRUE, help = "Whether to calculate CIs, default is TRUE")
  parser$add_argument('--aggregate', type = 'logical', default = FALSE, help = "Whether to aggregate age groups, default is FALSE")
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE, help = "Whether to exclude possible duplicate pairs, default is FALSE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
ci_flag <- args$ci
aggregate_flag <- args$aggregate
exclude_duplicates <- args$exclude_duplicates

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)

df_meta <- tbl(con,"metadata")
#Age aggregation if not done
if(aggregate_flag){
  local_meta <- df_meta %>%
    mutate(
      age_aggr = case_when(
        is.na(age_adj) ~ "NA",
        age_adj < 6 ~ "0-5y",
        age_adj < 12 ~ "6-11y",
        age_adj < 18 ~ "12-17y",
        age_adj < 26 ~ "18-25y",
        age_adj < 46 ~ "26-45y",
        age_adj < 66 ~ "46-65y",
        age_adj < 81 ~ "65-80y",
        .default = "81y+"
      ),
      school_cat = case_when(
        is.na(age_adj) ~ "NA",
        age_adj < 5 ~ "Pre-School",
        age_adj < 12 ~ "Primary School",
        age_adj < 18 ~ "Secondary School",
        age_adj < 65 ~ "Adult",
        .default = "Seniors"
      )
    ) %>% collect()
  DBI::dbWriteTable(con, "metadata", local_meta, overwrite = TRUE)
  df_meta <- tbl(con,"metadata")
}
MONTH_BUFFER = 1
# Make date sequences for the time series analysis
mid_start_date <- as.Date("2020-03-01")
mid_end_date <- as.Date("2024-10-01")
mid_date_vector <- seq(from = mid_start_date, to = mid_end_date,by = "4 weeks")
lb_date_vector <- mid_date_vector - 28 * MONTH_BUFFER
ub_date_vector <- mid_date_vector + 28 * MONTH_BUFFER

tic("Series anaylsis")
state_rr_series <- NULL
for(i in 1:length(lb_date_vector)){
  print(mid_date_vector[i])
  rr_series <- con %>%
    bind_pairs_exp("age_aggr",time_bounds = c(lb_date_vector[i],ub_date_vector[i]), exclude_duplicates = exclude_duplicates) %>%
    calculate_rr_matrix() %>%
    collect() %>%
    mutate(date = mean(c(lb_date_vector[i],ub_date_vector[i]))) %>%
    normalized_age_rr_diag() %>%
    normalized_age_rr_fixed(baseline_grp = "26-45y")

  if(i == 1){
    age_rr_series <- rr_series
  } else{
    age_rr_series <- bind_rows(age_rr_series,rr_series)
  }
}
toc()

fn_series <- paste0("results/",scenario,"/time_age/df_RR_by_time_age_series.tsv")
write_tsv(age_rr_series,file=fn_series)

# Country-specific age RR time series
tic("Age by country and time analysis")
age_country_time_rr <- NULL

# Get list of countries
countries <- df_meta %>%
  filter(!is.na(country)) %>%
  distinct(country) %>%
  pull(country)

# Create pairs with country information
pairs_country <- con %>%
  bind_pairs_exp("country", exclude_duplicates = exclude_duplicates) %>%
  rename(country.x = x) %>%
  rename(country.y = y) %>%
  mutate(sameCountry = (country.x == country.y))

for(i in 1:length(lb_date_vector)){
  print(paste0("Processing date: ", mid_date_vector[i]))

  # Bind pairs for the time period once
  pairs_time <- con %>%
    bind_pairs_exp("age_aggr", time_bounds = c(lb_date_vector[i], ub_date_vector[i]), exclude_duplicates = exclude_duplicates)
  pairs_time <- pairs_time %>%
    left_join(pairs_country, join_by(strain_1, strain_2))

  for(country in countries){
    # Filter pairs where both strains are from this country
    pair_subset <- pairs_time %>%
      filter(country.x == country & country.y == country) %>% collect()

    # Skip if no pairs for this country
    if(nrow(pair_subset) == 0) next

    # Calculate RR matrix for this country
    rr_country_time <- pair_subset %>%
      calculate_rr_matrix() %>%
      mutate(
        date = mean(c(lb_date_vector[i], ub_date_vector[i])),
        country = country
      ) %>%
      normalized_age_rr_diag(group_vars = c("date", "country")) %>%
      normalized_age_rr_fixed(baseline_grp = "26-45y", group_vars = c("date", "country"))

    # Bind to main dataframe
    if(is.null(age_country_time_rr)){
      age_country_time_rr <- rr_country_time
    } else {
      age_country_time_rr <- bind_rows(age_country_time_rr, rr_country_time)
    }
  }
}
toc()

fn_country_time <- paste0("results/", scenario, "/time_age/df_RR_by_time_age_countries.tsv")
write_tsv(age_country_time_rr, file = fn_country_time)


pairs_state <- con %>%
  bind_pairs_exp("division", exclude_duplicates = exclude_duplicates) %>%
  rename(division.x = x) %>%
  rename(division.y = y) %>%
  mutate(sameDiv = (division.x == division.y))

# School age analysis stratified by state and time
tic("School age by state and time analysis")
school_state_time_rr <- NULL

# Get list of states
states <- df_meta %>%
  filter(!is.na(division)) %>%
  distinct(division) %>%
  pull(division)

for(i in 1:length(lb_date_vector)){
  print(paste0("Processing date: ", mid_date_vector[i]))

  # Bind pairs for the time period once
  pairs_time <- con %>%
    bind_pairs_exp("school_cat", time_bounds = c(lb_date_vector[i], ub_date_vector[i]), exclude_duplicates = exclude_duplicates)
  pairs_time <- pairs_time %>%
    left_join(pairs_state,join_by(strain_1,strain_2)) 
  
  for(state in states){
    # Filter pairs where both strains are from this state
    pair_subset <- pairs_time %>%
      filter(division.x == state & division.y == state) %>% collect()

    # Skip if no pairs for this state
    if(nrow(pair_subset) == 0) next

    # Calculate RR matrix for this state
    rr_state_time <- pair_subset %>%
      calculate_rr_matrix() %>%
      mutate(
        date = mean(c(lb_date_vector[i], ub_date_vector[i])),
        state = state
      ) %>%
      normalized_age_rr_diag() %>%
      normalized_age_rr_fixed(baseline_grp = "Adult")

    # Bind to main dataframe
    if(is.null(school_state_time_rr)){
      school_state_time_rr <- rr_state_time
    } else {
      school_state_time_rr <- bind_rows(school_state_time_rr, rr_state_time)
    }
  }
}
toc()

fn_school_state_time <- paste0("results/", scenario, "/time_age/df_RR_by_school_state_time.tsv")
write_tsv(school_state_time_rr, file = fn_school_state_time)

# School age analysis stratified by state and academic year
tic("School age by state and academic year analysis")
school_state_ay_rr <- NULL

# Define academic years (Sept 1 to May 31)
academic_years <- data.frame(
  ay_label = c("2020-2021", "2021-2022", "2022-2023", "2023-2024"),
  start_date = as.Date(c("2020-09-01", "2021-09-01", "2022-09-01", "2023-09-01")),
  end_date = as.Date(c("2021-05-31", "2022-05-31", "2023-05-31", "2024-05-31"))
)

for(ay_idx in 1:nrow(academic_years)){
  ay_label <- academic_years$ay_label[ay_idx]
  start_date <- academic_years$start_date[ay_idx]
  end_date <- academic_years$end_date[ay_idx]

  print(paste0("Processing academic year: ", ay_label))

  # Bind pairs for the time period once
  pairs_time <- con %>%
    bind_pairs_exp("school_cat", time_bounds = c(start_date, end_date), exclude_duplicates = exclude_duplicates)
  pairs_time <- pairs_time %>%
    left_join(pairs_state,join_by(strain_1,strain_2)) 

  for(state in states){
    # Filter pairs where both strains are from this state
    pairs_subset <- pairs_time %>%
      filter(division.x == state & division.y == state) %>%
      collect()

    # Skip if no pairs for this state
    if(nrow(pairs_subset) == 0) next

    # Calculate RR matrix for this state
    rr_state_ay <- pairs_subset %>%
      calculate_rr_matrix() %>%
      mutate(
        academic_year = ay_label,
        state = state
      ) %>%
      normalized_age_rr_diag() %>%
      normalized_age_rr_fixed(baseline_grp = "Adult")

    # Bind to main dataframe
    if(is.null(school_state_ay_rr)){
      school_state_ay_rr <- rr_state_ay
    } else {
      school_state_ay_rr <- bind_rows(school_state_ay_rr, rr_state_ay)
    }
  }
}
toc()

fn_school_state_ay <- paste0("results/", scenario, "/time_age/df_RR_by_school_state_ay.tsv")
write_tsv(school_state_ay_rr, file = fn_school_state_ay)

DBI::dbDisconnect(con, shutdown = TRUE)
