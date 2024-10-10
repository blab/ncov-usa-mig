#File: age_state_analysis.R
#Author(s): Amin Bemanian
#Date: 8/22/24
#Description: Calculate the age relative risks within each state and stratified by adjacency
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--threads: Number of threads for the parallelization, default is 1 if no parallelization
#Notes:
#Requires running age_analysis state_analysis first, due to conserving the same dataframe for pairs

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(doParallel)
library(foreach)

source('./scripts/bind_pairs_exp.R')
source("./scripts/calculate_rr_matrix.R")

STATE_DISTANCES <- fread("../data/nb_dist_states.tsv")


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "USA", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--threads', type = 'integer', default = 1, help = "How many threads for parallelization")
  return(parser$parse_args())
}

args <- collect_args()

scenario <- args$scenario
threads <- args$threads
N_SAMP <- 20 #Ideally do 1000 in final

registerDoParallel(cores = threads)


#x and y refer to state/division in pairs
get_adj_status <- function(x,y){
  d_nb <- STATE_DISTANCES[state_x == x & state_y == y, nb_dist]
  case_when(
    d_nb == 0 ~ "Within",
    d_nb == 1 ~ "Adjacent",
    TRUE ~ "Non-adjacent"
  )
}

tryCatch(
  {
    df_pairs_with_age <- fread(paste0("results/",scenario,"/df_pairs_with_age.tsv"))
    df_pairs_with_state  <- fread(paste0("results/",scenario,"/df_pairs_with_state.tsv"))
    df_pairs_age_state <- left_join(df_pairs_with_age,df_pairs_with_state)
  },
  error = function(e){
    print("Missing pair datasets joined with age and state. Run those analyses prior.")
    stop()
  }
)

#Assign adjacency status
df_pairs_age_state$adj_status <- ""
LIST_STATES <- union(unique(df_pairs_age_state$X_division),unique(df_pairs_age_state$Y_division))
#Far quicker instead of iterating over each pair and running get_adj_status
for(x in LIST_STATES){
  for(y in LIST_STATES){
    adj <- get_adj_status(x,y)
    df_pairs_age_state[X_division == x & Y_division == y]$adj_status <- adj
  }
}

#First make an analysis within each state for Within state pairs
df_age_RRs_by_state <- tibble()

for (state in LIST_STATES) {
  df_p_age <- df_pairs_age_state[X_division == state & adj_status == 'Within',]
  tryCatch(
    {
      ages_rr <- calculate_rr_matrix(df_p_age,"age_class",return_ci = TRUE,n_subsample=N_SAMP) %>%
        mutate(state_name = state)
      df_age_RRs_by_state <- rbind(df_age_RRs_by_state,ages_rr)
    },
    error = function(e){
      print(paste0("Failed to do age RR analysis on: ",state))
    }
  )
}

fn_age_RRs_by_state <- paste0("results/",scenario,"/df_RR_by_age_state.tsv")
readr::write_tsv(df_age_RRs_by_state,file = fn_age_RRs_by_state)


df_age_RRs_by_adj <- tibble()
for(status in c("Within","Adjacent","Non-adjacent")){
  ages_rr <- df_pairs_age_state %>% 
    filter(adj_status == status) %>%
    calculate_rr_matrix(exp_var="age_class",return_ci = TRUE) %>%
    mutate(adj_status = status)
  df_age_RRs_by_adj <- rbind(df_age_RRs_by_adj,ages_rr)
}

fn_age_RRs_by_adj <- paste0("results/",scenario,"/df_RR_by_age_adj.tsv")
readr::write_tsv(df_age_RRs_by_adj,file = fn_age_RRs_by_adj)