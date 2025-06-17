#File: census_div_analysis.R
#Author(s): Amin Bemanian
#Date: 8/19/24
#Description: Calculate the relative risks across states
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--threads: Number of threads for the parallelization, default is 1 if no parallelization

library(tidyverse)
library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(duckdb)
library(dbplyr)


source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")

STATE_DISTANCES <- fread("./data/nb_dist_states.tsv")

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

div_rr <- con %>% 
  bind_pairs_exp("bea_reg") %>%
  calculate_rr_matrix() %>%
  collect()

if(ci_flag){
  div_rr_ci <- calculate_rr_ci(con,"bea_reg")
  div_rr <- inner_join(div_rr,div_rr_ci,by=join_by(x,y))
}

fn_rr <- paste("results/",scenario,"/df_RR_by_census_div.tsv",sep="")
readr::write_tsv(div_rr,file=fn_rr)
DBI::dbDisconnect(con,shutdown=TRUE)
print("Successfully finished census division analysis!")
