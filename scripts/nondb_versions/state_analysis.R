#File: state_analysis.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Calculate the relative risks across states
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--threads: Number of threads for the parallelization, default is 1 if no parallelization

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)

source('./scripts/bind_pairs_exp.R')
source("./scripts/calculate_rr_matrix.R")

STATE_DISTANCES <- fread("./data/nb_dist_states.tsv")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--threads', type = 'integer', default = 1, help = "How many threads for parallelization")
  parser$add_argument('--ci', type = 'logical', default = TRUE, help = "Whether to calculate CIs, default is TRUE")
  return(parser$parse_args())
}

args <- collect_args()

scenario <- args$scenario
threads <- args$threads
ci_flag <- args$ci

registerDoParallel(cores = threads)

fn_df_meta <- paste0("./results/",scenario,"/df_meta_clust.tsv")
fn_pairs <- paste0('./results/',scenario,"/df_pairs_clean.tsv") 

fn_pairs_with_state <- paste0("results/",scenario,"/df_pairs_with_state.tsv")

Sys.time()
df_pairs_with_state <- bind_pairs_exp(fn_m = fn_df_meta,
               fn_p = fn_pairs,
               exp_var = "division",
               fn_o = fn_pairs_with_state)

state_rr <- calculate_rr_matrix(df_pairs_with_state,"division",return_ci = ci_flag)
state_rr <- state_rr %>% left_join( #Attach the adjacency and distance variables
  STATE_DISTANCES,
  by=c(
    "x"="state_x",
    "y"="state_y"
  )
)
Sys.time()

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
readr::write_tsv(state_rr,file=fn_rr)

print("Successfully finished state analysis!")
