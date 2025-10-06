library(tidyverse)
library(dbplyr)
library(duckdb)
library(tictoc)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")

collect_args <- function(){
  parser <- argparse::ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--ci', type = 'logical', default = TRUE, help = "Whether to calculate CIs, default is TRUE")
  parser$add_argument('--aggregate', type = 'logical', default = FALSE, help = "Whether to aggregate age groups, default is FALSE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
ci_flag <- args$ci
aggregate_flag <- args$aggregate

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
        age_adj < 18 ~ "6-17y",
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
MONTH_BUFFER = 3 
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
  if(i == 1){
    age_rr_series <- con %>% 
      bind_pairs_exp("age_aggr",time_bounds = c(lb_date_vector[i],ub_date_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(lb_date_vector[i],ub_date_vector[i])))
  } else{
    rr_series <- con %>% 
      bind_pairs_exp("age_aggr",time_bounds = c(lb_date_vector[i],ub_date_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(lb_date_vector[i],ub_date_vector[i])))
    age_rr_series <- bind_rows(age_rr_series,rr_series)
  } 
}
toc()

fn_series <- paste0("results/",scenario,"/df_RR_by_age_time_series.tsv")
write_tsv(age_rr_series,file=fn_series)

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
    bind_pairs_exp("school_cat", time_bounds = c(lb_date_vector[i], ub_date_vector[i])) %>%
    collect()

  for(state in states){
    # Filter pairs where both strains are from this state
    pairs_state <- pairs_time %>%
      filter(division_x == state & division_y == state)

    # Skip if no pairs for this state
    if(nrow(pairs_state) == 0) next

    # Calculate RR matrix for this state
    rr_state_time <- pairs_state %>%
      calculate_rr_matrix() %>%
      mutate(
        date = mean(c(lb_date_vector[i], ub_date_vector[i])),
        state = state
      )

    # Bind to main dataframe
    if(is.null(school_state_time_rr)){
      school_state_time_rr <- rr_state_time
    } else {
      school_state_time_rr <- bind_rows(school_state_time_rr, rr_state_time)
    }
  }
}
toc()

fn_school_state_time <- paste0("results/", scenario, "/df_RR_by_school_state_time.tsv")
write_tsv(school_state_time_rr, file = fn_school_state_time)

DBI::dbDisconnect(con, shutdown = TRUE)
