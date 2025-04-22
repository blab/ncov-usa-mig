#File: clean_data.R
#Author(s): Amin Bemanian
#Date: 8/8/24
#Description: Simple cleaning of dataset to make a curated dataframe that combines the cluster data with the cluster allocations

library(argparse)
library(dplyr)
library(data.table)
library(dbplyr)
library(duckdb)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Dataframe with pairs of identical sequences')
  return(parser$parse_args())
}

#To clean age function
clean_age <- function(a){
  UB <- 125 #Maximum age, oldest documented humans in history are between 117-122 in age
  a <- suppressWarnings(as.numeric(a))  # Convert safely, suppress warnings
  a <- ifelse(!is.na(a) & a >= 0 & a < UB, as.integer(a), NA_integer_)  # Ensure integer output
  return(a)
}

#Aggregate the ages for the RR analysis, note that age should be cleaned prior
aggregate_age <- function(a, BIN_SIZE = 5){
  if (scenario == 'USA') {  #For USA will use individual age bins
    return(case_when(
      is.na(a) ~ 'NA',  # Handle missing values
      a < 99 ~ paste0(sprintf("%02d", a), "y"),  # Format ages 0-98 as "XXy"
      TRUE ~ '99y+'  # Group ages 99+
    ))
  } else { #For all other scenarios will likely need larger bins, default of 5
    return(case_when(
      is.na(a) ~ 'NA',
      a < (90 %/% BIN_SIZE) * BIN_SIZE ~ paste0(
        sprintf("%02d", (a %/% BIN_SIZE) * BIN_SIZE), "-",
        sprintf("%02d", (a %/% BIN_SIZE) * BIN_SIZE + (BIN_SIZE - 1)), "y"
      ), # Group into BIN_SIZE-year bins
      TRUE ~ '90y+'  # Group ages 90+
    ))
  }
}


args <- collect_args()
scenario <- args$scenario

dir.create(paste0("results/",scenario)) #Make a directory in case it doesn't already exist

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)
df_meta <- tbl(con,"metadata")
total_seq_n <- df_meta |> summarize(n()) |> collect()
print(paste0("Total number of sequences initially: ",total_seq_n))

df_meta <- df_meta %>%
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
  filter(census_reg %in% c("West","Midwest","South","Northeast","Hawaii","Alaska")) %>% #Limits to 50 States + US
  filter(!is.na(division))  %>%
  collect() %>%
  mutate(age_adj = clean_age(age)) %>%
  mutate(age_adj = as.integer(age_adj)) %>% # Explicitly cast to integer
  mutate(age_class = aggregate_age(age_adj))


dbWriteTable(conn = con,
        value = df_meta,
        name = "metadata",
        overwrite=TRUE
)

filtered_seq_n <- df_meta |> summarize(n()) 
print(paste0("Number of sequences post-filtering: ",filtered_seq_n))


#Prune pairs to only use strains that were included in the final metadata list
#Manually writing SQL for this step as the dbplyr grammar struggles with cross table comparisons
# df_pairs <- sql("
# SELECT * FROM pairs
#   WHERE strain_1 IN (SELECT strain FROM metadata)
#   AND strain_2 IN (SELECT strain FROM metadata)
# ") %>% dbGetQuery(con,.)
# print("Successfully got pairs")
#   
# dbWriteTable(conn = con,
#                value = df_pairs,
#                name = "pairs",
#                overwrite=TRUE
# )
    
DBI::dbDisconnect(con, shutdown=TRUE)
