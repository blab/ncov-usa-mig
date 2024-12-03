#File: age_analysis.R
#Author(s): Amin Bemanian
#Date: 8/15/24
#Description: Simple cleaning of dataset to make a curated dataframe that combines the cluster data with the cluster allocations
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--threads: Number of threads for the parallelization, default is 1 if no parallelization

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(duckdb)
library(dbplyr)

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

ages_rr <- con %>% 
  bind_pairs_exp("age_class") %>%
  calculate_rr_matrix() %>%
  collect()

if(ci_flag){
  ages_rr_ci <- calculate_rr_ci(con,"age_class")
  ages_rr <- inner_join(ages_rr,ages_rr_ci,by=join_by(x,y))
}

fn_rr <- paste("results/",scenario,"/df_RR_by_age_class.tsv",sep="")
readr::write_tsv(ages_rr,file=fn_rr)
DBI::dbDisconnect(con)
print("Successfully finished age analysis!")