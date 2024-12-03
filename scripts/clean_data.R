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
  a <- as.numeric(a) #Force age to be a number
  return(ifelse(!is.na(a) & a >= 0 & a < UB,a,NA))
}

#Aggregate the ages for the RR analysis, note that age should be cleaned prior
aggregate_age <- function(a){
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
                     a <= 9 ~ '05-09y',
                     a <= 14 ~ '10-14y',
                     a <= 19 ~ '15-19y',
                     a <= 24 ~ '20-24y',
                     a <= 29 ~ '25-29y',
                     a <= 34 ~ '30-34y',
                     a <= 39 ~ '34-39y',
                     a <= 44 ~ '40-44y',
                     a <= 49 ~ '45-49y',
                     a <= 54 ~ '50-54y',
                     a <= 59 ~ '55-59y',
                     a <= 64 ~ '60-64y',
                     a <= 69 ~ '65-69y',
                     a <= 74 ~ '70-74y',
                     a <= 79 ~ '75-79y',
                     a <= 84 ~ '80-84y',
                     a <= 89 ~ '85-89y',
                     TRUE ~ '90y+'))
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
df_pairs <- sql("
SELECT * FROM pairs
  WHERE strain_1 IN (SELECT strain FROM metadata)
  AND strain_2 IN (SELECT strain FROM metadata)
") %>% dbGetQuery(con,.)

  
dbWriteTable(conn = con,
               value = df_pairs,
               name = "pairs",
               overwrite=TRUE
)
    
DBI::dbDisconnect(con, shutdown=TRUE)
