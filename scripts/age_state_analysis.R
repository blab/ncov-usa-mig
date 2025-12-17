#File: age_state_analysis.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Calculate the age relative risks within each state and stratified by adjacency
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--threads: Number of threads for the parallelization, default is 1 if no parallelization
#Notes:
#Requires running age_analysis state_analysis first, due to conserving the same dataframe for pairs

library(tidyverse)
library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(duckdb)
library(dbplyr)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")

STATE_DISTANCES <- fread("data/nb_dist_states.tsv")
REGION_LIST <- fread("data/us_states_regions.csv")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--ci', type = 'logical', default = FALSE, help = "Whether to calculate CIs, default is TRUE")
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE, help = "Whether to exclude possible duplicate pairs, default is FALSE")
  return(parser$parse_args())
}


args <- collect_args()
ci_flag <- args$ci
scenario <- args$scenario
exclude_duplicates <- args$exclude_duplicates

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)

pairs_state <- con %>%
  bind_pairs_exp("division", exclude_duplicates = exclude_duplicates) %>%
  rename(state.x = x) %>%
  rename(state.y = y) %>%
  mutate(sameState = (state.x == state.y)) 

pairs_region <- con %>%
  bind_pairs_exp("bea_reg", exclude_duplicates = exclude_duplicates) %>%
  rename(reg.x = x) %>%
  rename(reg.y = y) %>%
  mutate(sameRegion = (reg.x == reg.y)) 

pairs_state_reg <- inner_join(pairs_state,pairs_region)

pairs_age <- con %>%
  bind_pairs_exp("age_class", exclude_duplicates = exclude_duplicates) 

pairs_state_age <- inner_join(pairs_state_reg,pairs_age)

#Need for CI subsampling
meta_age <- con %>%
  tbl("metadata") %>%
  filter(!is.na(age_class)) %>%
  filter(age_class != 'NA')

#This analysis needs a unique CI calculation since it's being stratified over several add'l variables(state, sameState, sameRegion)
#pair_list - full list, will be pairs_state_age
#meta_tbl - Metadata table to subset sequences from
#specific_state - name of a specific state to stratify
#TODO 12/11 - IMPLEMENT STRATIFIED CI FOR REGIONS
calculate_rr_ci <- function(con,meta_tbl,specific_state=NULL,samp_cov=0.8,k=200,interval_width=0.95,exclude_duplicates=FALSE){
  RR_dist_same <- matrix()
  RR_dist_diff <- matrix()
  for(x in 1:k){
    pairs_state_samp <- con %>%
      bind_pairs_exp("division", exclude_duplicates = exclude_duplicates) %>%
      rename(state.x = x) %>%
      rename(state.y = y) %>%
      mutate(sameState = (state.x == state.y))
    if(!is.null(specific_state)){
      pairs_state_samp <- pairs_state_samp %>% filter(state.x == specific_state | state.y == specific_state)
    }
    pairs_age_samp <- con %>%
      bind_pairs_exp("age_class",sub_samp = TRUE,samp_cov, exclude_duplicates = exclude_duplicates)
    pairs_combo_samp <- inner_join(pairs_state_samp,pairs_age_samp)
    
    RR_samp_same <- pairs_combo_samp %>%
      filter(sameState) %>%
      calculate_rr_matrix() %>%
      collect()
    RR_samp_diff <- pairs_combo_samp %>%
      filter(!sameState) %>%
      calculate_rr_matrix() %>%
      collect()
    
    if(x == 1){
      CI_keys <- RR_samp_same %>% select(c(x,y))
      RR_dist_same <- RR_samp_same$RR
      RR_dist_diff <- RR_samp_diff$RR
    }else{
      RR_dist_same <- cbind(RR_dist_same,RR_samp_same$RR)
      RR_dist_diff <- cbind(RR_dist_diff,RR_samp_diff$RR)
    }
  }
  ci_lb <- apply(RR_dist_same,1,quantile,probs=(1-interval_width)/2)
  ci_ub <- apply(RR_dist_same,1,quantile,probs=(1+interval_width)/2)
  CI_same <- cbind(CI_keys,ci_lb,ci_ub)
  ci_lb <- apply(RR_dist_diff,1,quantile,probs=(1-interval_width)/2)
  ci_ub <- apply(RR_dist_diff,1,quantile,probs=(1+interval_width)/2)
  CI_diff <- cbind(CI_keys,ci_lb,ci_ub)

  return(list(CI_same,CI_diff))
}


ages_rr_same_state <- pairs_state_age %>%
  filter(sameState) %>%
  calculate_rr_matrix() %>%
  collect() %>%
  mutate(sameState=TRUE, sameRegion=TRUE)

ages_rr_same_region <- pairs_state_age %>%
  filter(!sameState && sameRegion) %>%
  calculate_rr_matrix() %>%
  collect() %>%
  mutate(sameState = FALSE, sameRegion=TRUE)

ages_rr_diff_region <- pairs_state_age %>%
  filter(!sameRegion) %>%
  calculate_rr_matrix() %>%
  collect() %>%
  mutate(sameState=FALSE, sameRegion=FALSE) 

#TODO Implement CI for same_region
if(ci_flag){
  ci_list <- calculate_rr_ci(con,meta_age,samp_cov=0.8,exclude_duplicates=exclude_duplicates)
  ages_rr_same <- inner_join(ages_rr_same,ci_list[[1]],by=join_by(x,y))
  ages_rr_diff <- inner_join(ages_rr_diff,ci_list[[2]],by=join_by(x,y))
}

#Combine different and same state RRs
age_state_rr <- rbind(ages_rr_same_state,
  ages_rr_same_region,
  ages_rr_diff_region) 
fn_rr <- paste("results/",scenario,"/df_RR_by_age_state.tsv",sep="")
readr::write_tsv(age_state_rr,file=fn_rr)
DBI::dbDisconnect(con,shutdown=TRUE)
print("Successfully finished age/state stratified analysis!")

