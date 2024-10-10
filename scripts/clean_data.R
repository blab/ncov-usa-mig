#File: clean_data.R
#Author(s): Amin Bemanian
#Date: 8/8/24
#Description: Simple cleaning of dataset to make a curated dataframe that combines the cluster data with the cluster allocations

library(argparse)
library(dplyr)
library(data.table)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Dataframe with pairs of identical sequences')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

dir.create(paste0("results/",scenario)) #Make a directory in case it doesn't already exist

fn_meta <- paste0("./data/",scenario,"/metadata/metadata_",scenario,".tsv")
fn_cluster <- paste0("./data/",scenario,"/distance_aggregated/combined_df_cluster_alloc_",scenario,".tsv")
fn_pairs <- paste0('./data/',scenario,"/distance_aggregated/combined_df_identical_pairs_",scenario,".tsv") 

df_meta <- fread(fn_meta)
print(paste0("Total number of sequences initially:",nrow(df_meta)))

if(scenario=="USA"){ #For the full USA scenario, cluster allocation code is not working, will bypass the allocation since we are not using for now
  df_cluster_full <- df_meta #Just assign it to be the metadata without the clusters
} else{ #All other scenarios, goal is to have these steps included into the full set at a later date
  df_cluster_alloc <- fread(fn_cluster)

  #Make a unique cluster ID to simplify analysis
  df_cluster_alloc <- df_cluster_alloc %>% mutate(cluster_pango_id = paste(Nextclade_pango,cluster_id,sep = "_")) %>%
    select(strain,cluster_pango_id)

  #Combine the tables and do some other filtering/cleaning of the metadata
  df_cluster_full <- left_join(df_meta,df_cluster_alloc,by = "strain") 
}


df_cluster_full <- df_cluster_full %>%
  mutate(strain_1 = strain) %>% #Need to create dummy columns for tsv-utils join to the pairs in bind_pairs_exp
  mutate(strain_2 = strain) %>%
  filter(host=="Human") %>% #Limit to human cases only
  filter(coverage >= 0.9) %>% #90% or higher coverage
  mutate(division = case_when(
    division %in% c("Washington DC") ~ "District of Columbia",#To make it consistent with the Census Bureau convention
    TRUE ~ division 
    )
  ) %>%
  mutate(census_div = case_when(  #Attach census divisions to metadata
    division %in% c("Washington","Oregon","California") ~ "Pacific",
    division %in% c("Nevada","Idaho","Montana","Wyoming","Utah","Colorado","Arizona","New Mexico") ~ "Mountain",
    division %in% c("North Dakota","South Dakota","Nebraska","Kansas","Minnesota","Iowa","Missouri") ~ "West North Central",
    division %in% c("Wisconsin","Illinois","Indiana","Michigan","Ohio") ~ "East North Central",
    division %in% c("Oklahoma","Arkansas","Texas","Louisiana") ~ "West South Central",
    division %in% c("Kentucky","Tennessee","Mississippi","Alabama") ~ "East South Central",
    division %in% c("Florida","Georgia","South Carolina","North Carolina","Virginia","West Virginia","District of Columbia","Maryland","Delaware") ~ "South Atlantic",
    division %in% c("New York","Pennsylvania","New Jersey") ~ "Middle Atlantic",
    division %in% c("Connecticut","Rhode Island","Massachusetts","Vermont","New Hampshire","Maine") ~ "New England",
    division %in% c("Alaska") ~ "Alaska",
    division %in% c("Hawaii") ~ "Hawaii",
    division %in% c("Guam","American Samoa") ~ "Pacific Territories",
    division %in% c("Puerto Rico","Virgin Islands") ~ "Caribbean Territories",
    TRUE ~ "Invalid"
    )
  ) %>%
  mutate(census_reg = case_when( #Aggregate to census regions
    census_div %in% c("Pacific","Mountain") ~ "West",
    census_div %in% c("West North Central","East North Central") ~ "Midwest",
    census_div %in% c("West South Central","East South Central","South Atlantic") ~ "South",
    census_div %in% c("Middle Atlantic","New England") ~ "Northeast",
    TRUE ~ census_div
    )
  ) %>%
  filter(census_reg %in% c("West","Midwest","South","Northeast","Hawaii","Alaska")) %>%
  filter(!is.na(division))

print(paste0("Number of sequences post-filtering:",nrow(df_cluster_full)))

fn_output <- paste0("./results/",scenario,"/df_meta_clust.tsv")
readr::write_tsv(df_cluster_full,fn_output)

STRAIN_LIST <- df_cluster_full$strain 
df_pairs <- fread(fn_pairs)
print(paste0("Number of pairs pre-filtering: ",nrow(df_pairs)))

#Filter the pairs to only include pairs where both sequences are inclcuded in the metadata strain list
df_pairs_clean <- df_pairs %>% filter(strain_1 %in% STRAIN_LIST & strain_2 %in% STRAIN_LIST)
print(paste0("Number of pairs post-filtering: ",nrow(df_pairs_clean)))

fn_pairs_out <- paste0("./results/",scenario,"/df_pairs_clean.tsv")
readr::write_tsv(df_pairs_clean,fn_pairs_out)