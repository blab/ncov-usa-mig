#File: census_div_analysis.R
#Author(s): Amin Bemanian
#Date: 8/19/24
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

df_meta <- fread(fn_df_meta) 
df_pairs <- fread(fn_pairs)

fn_pairs_with_census_div <- paste0("results/",scenario,"/df_pairs_with_census_div.tsv")

df_pairs_with_census_div <- bind_pairs_exp(fn_m = fn_df_meta,
               fn_p = fn_pairs,
               exp_var = "census_div",
               fn_o = fn_pairs_with_census_div)

#Identical pairs analysis
div_rr <- calculate_rr_matrix(df_pairs_with_census_div,"census_div",return_ci = ci_flag)
fn_rr <- paste("results/",scenario,"/df_RR_by_census_div.tsv",sep="")
readr::write_tsv(div_rr,file=fn_rr)

