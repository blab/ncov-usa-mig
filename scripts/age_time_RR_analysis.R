library(tidyverse)
library(dbplyr)
library(duckdb)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--ci', type = 'logical', default = TRUE, help = "Whether to calculate CIs, default is TRUE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
ci_flag <- args$ci

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)

df_meta <- tbl(con,"metadata")
if(!("age_aggr" %in% colnames(df_meta))){
  df_meta <- df_meta %>%
    mutate(age_aggr = case_when(
        age < 6 ~ "0-5y",
        age < 18 ~ "6-17y",
        age < 26 ~ "18-25y",
        age < 46 ~ "26-45y",
        age < 66 ~ "65+y"
    ))
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
  if(i == 1){
    state_rr_series <- con %>% 
      bind_pairs_exp("division",time_bounds = c(lb_date_vector[i],ub_date_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(lb_date_vector[i],ub_date_vector[i])))
  } else{
    rr_series <- con %>% 
      bind_pairs_exp("division",time_bounds = c(lb_date_vector[i],ub_date_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(lb_date_vector[i],ub_date_vector[i])))
    state_rr_series <- bind_rows(state_rr_series,rr_series)
  } 
}
toc()