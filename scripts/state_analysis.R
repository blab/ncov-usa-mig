#File: state_analysis.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Calculate the relative risks across states
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--threads: Number of threads for the parallelization, default is 1 if no parallelization

library(argparse)
library(tidyverse)
library(data.table)
library(ggplot2)
library(duckdb)
library(dbplyr)


source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")

STATE_DISTANCES <- fread("data/nb_dist_states.tsv")
CBSA_DISTANCE <- fread("data/state_cbsa_dist.tsv")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--ci', type = 'logical', default = TRUE, help = "Whether to calculate CIs, default is TRUE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
ci_flag <- args$ci

#Debugging shortcuts
#scenario <- "CAM_10"
#ci_flag <- FALSE

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)

state_rr <- con %>% 
  bind_pairs_exp("division") %>%
  calculate_rr_matrix() %>%
  collect()

if(ci_flag){
  state_rr_ci <- calculate_rr_ci(con,"division")
  state_rr <- inner_join(state_rr,state_rr_ci,by=join_by(x,y))
}

state_rr <- state_rr %>% left_join( #Attach the adjacency and distance variables
  STATE_DISTANCES,
  by=c(
    "x"="state_x",
    "y"="state_y"
  )
) %>% 
  left_join(CBSA_DISTANCE,by=join_by(x,y))

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
print(fn_rr)
readr::write_tsv(state_rr,file=fn_rr)
DBI::dbDisconnect(con,shutdown=TRUE)
print("Successfully finished state analysis!")
