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
library(doParallel)
library(foreach)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")


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


df_meta_clust <- fread(fn_df_meta)
df_pairs <- fread(fn_pairs)

#Make sure the age is within a human feasible range
clean_age <- function(a){
  UB <- 125 #Maximum age, oldest documented humans in history are between 117-122 in age
  a <- as.numeric(a) #Force age to be a number
  return(ifelse(!is.na(a) & a >= 0 & a < UB,a,NA))
}

#Aggregate the ages for the RR analysis, note that age should be cleaned prior
aggregate_age <- function(a){
  # return(case_when(is.na(as.numeric(a)) ~ 'NA', #Just for safety checking
  #           a <= 9 ~ '00-09y', 
  #           a <= 19 ~ '10-19y',
  #           a <= 29 ~ '20-29y', 
  #           a <= 39 ~ '30-39y',
  #           a <= 49 ~ '40-49y', 
  #           a <= 59 ~ '50-59y',
  #           a <= 69 ~ '60-69y',
  #           a <= 79 ~ '70-79y',
  #           TRUE ~ '80y+'))
  
  # Alternate age structure
  if(scenario == 'USA'){ #Divide by every 3 years
    return(case_when(is.na(as.numeric(a)) ~ 'NA', #Just for safety checking
            a <= 2 ~ '00-02y', 
            a <= 5 ~ '03-05y',
            a <= 8 ~ '06-08y', 
            a <= 11 ~ '09-11y',
            a <= 14 ~ '12-14y', 
            a <= 17 ~ '15-17y',
            a <= 20 ~ '18-20y',
            a <= 23 ~ '21-23y',
            a <= 26 ~ '24-26y',
            a <= 29 ~ '27-29y',
            a <= 32 ~ '30-32y',
            a <= 35 ~ '33-35y',
            a <= 38 ~ '36-38y',
            a <= 41 ~ '39-41y',
            a <= 44 ~ '42-44y',
            a <= 47 ~ '45-47y',
            a <= 50 ~ '48-50y',
            a <= 53 ~ '51-53y',
            a <= 56 ~ '54-56y',
            a <= 59 ~ '57-59y',
            a <= 62 ~ '60-62y',
            a <= 65 ~ '63-65y',
            a <= 68 ~ '66-68y',
            a <= 71 ~ '69-71y',
            a <= 74 ~ '72-74y',
            a <= 77 ~ '75-77y',
            a <= 80 ~ '78-80y',
            a <= 83 ~ '80-83y',
            a <= 86 ~ '84-86y',
            a <= 89 ~ '86-89y',
            TRUE ~ '90y+'))
  }else{
    return(case_when(is.na(as.numeric(a)) ~ 'NA', #Just for safety checking
            a <= 4 ~ '00-04y', 
            a <= 11 ~ '05-11y',
            a <= 19 ~ '12-19y', 
            a <= 34 ~ '20-34y',
            a <= 49 ~ '35-49y', 
            a <= 64 ~ '50-64y',
            a <= 79 ~ '65-79y',
            TRUE ~ '80y+'))
  }
  
}


#Cleaning and aggregation of ages in data set
df_meta_clust <- df_meta_clust %>%
  mutate(age_adj = clean_age(age)) %>%
  mutate(age_class = aggregate_age(age_adj))

perc_age_valid <- 100 * sum(!is.na(df_meta_clust$age_adj))/nrow(df_meta_clust) 
print(
  paste("The percent of sequences with a valid age is ",round(perc_age_valid,digits=1),"%",sep="")
)
print("Table of ages (NA includes missing and invalid ages):")
print(table(df_meta_clust$age_class))

#Filter the datasets to only include valid ages, otherwise affects the RR calculation
df_meta_clust <- df_meta_clust %>%
  filter(!is.na(age_adj))
 
valid_strains <- df_meta_clust %>% pull(strain)
df_pairs_valid <- df_pairs %>% 
  filter(strain_1 %in% valid_strains & strain_2 %in% valid_strains)
perc_pairs_valid <- 100 * nrow(df_pairs_valid)/nrow(df_pairs)

print(
  paste("The percent of pairs with both sequences having valid ages is ",round(perc_pairs_valid,digits=1),"%",sep="")
)

rm(df_pairs)
gc()

fn_df_meta_age <- paste0("results/",scenario,"/metadata_age_clean.tsv")
fn_df_pair_clean <- paste0("results/",scenario,"/df_pairs_age_clean.tsv") #This file is just for the pairs with invalid ages removed

readr::write_tsv(df_meta_clust,file=fn_df_meta_age)
readr::write_tsv(df_pairs_valid,file=fn_df_pair_clean)

#Do some cleaning before entering the RR calculation step
rm(df_pairs_valid)
rm(df_meta_clust)
rm(valid_strains)
gc() #Helpful since R will be calling on tsv-utils to do the binding and it can't clean up R memory

#Calculation of RR by identical pairs

fn_pairs_with_age <- paste0("results/",scenario,"/df_pairs_with_age.tsv")
df_pairs_with_age <- bind_pairs_exp(fn_m=fn_df_meta_age,
                                    fn_p=fn_df_pair_clean,
                                    exp_var="age_class",
                                    fn_o=fn_pairs_with_age)

ages_rr <- calculate_rr_matrix(df_pairs_with_age,"age_class",return_ci = ci_flag)


fn_rr <- paste("results/",scenario,"/df_RR_by_age_class.tsv",sep="")
readr::write_tsv(ages_rr,file=fn_rr)

print("Successfully finished age analysis!")
