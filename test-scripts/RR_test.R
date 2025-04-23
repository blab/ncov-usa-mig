#File: RR_test.R
#Author(s): Amin Bemanian
#Date: 4/20/24
#Description: Creates a test database for debugging the RR calculaltion code

library(tidyverse)
library(data.table)
library(dbplyr)
library(duckdb)
library(ape)

source("scripts/calculate_rr_matrix.R")

FLAG_MAKE_PAIRS <- FALSE #If true will make DB on fly, otherwise takes in fasta + metadata
FASTA_FILE <- "test_data/synthetic_data/synthetic-fasta2.fasta"
META_FILE <- "test_data/synthetic_data/synthetic-metadata.csv"

#Make the test table
if(FLAG_MAKE_PAIRS){
  df_mat <- data.frame(
    x = character(),
    y = character(),
    n_xy = numeric()
  )
  
  df_mat <- df_mat %>%
    bind_rows(tibble(x = 'A', y = 'A', n_xy = 500000)) %>%
    bind_rows(tibble(x = 'A', y = 'B', n_xy = 5)) %>%
    bind_rows(tibble(x = 'A', y = 'C', n_xy = 5)) %>%
    bind_rows(tibble(x = 'B', y = 'A', n_xy = 5)) %>%
    bind_rows(tibble(x = 'B', y = 'B', n_xy = 500)) %>%
    bind_rows(tibble(x = 'B', y = 'C', n_xy = 5)) %>%
    bind_rows(tibble(x = 'C', y = 'A', n_xy = 5)) %>%
    bind_rows(tibble(x = 'C', y = 'B', n_xy = 5)) %>%
    bind_rows(tibble(x = 'C', y = 'C', n_xy = 50))
  
  test_pairs <- df_mat %>% uncount(n_xy)
} else{
  print(paste('Reading FASTA file:', FASTA_FILE))
  sequence_data <- read.FASTA(FASTA_FILE)
  print(paste('Reading metadata file:', META_FILE))
  metadata_data <- read.csv(META_FILE) %>% as_tibble()
  
  dist_mat <- dist.dna(sequence_data, model = 'N', as.matrix = T, pairwise.deletion = T)
  test_pairs <- dist_mat %>% 
    as_tibble() %>% 
    mutate(label_tip_1 = rownames(dist_mat)) %>% 
    pivot_longer(cols = -'label_tip_1', names_to = 'label_tip_2', values_to = 'n_mutations') %>% 
    filter(label_tip_1 != label_tip_2) %>%
    filter(n_mutations == 0)
  
  test_pairs <- test_pairs %>%
    left_join(metadata_data,
              by=join_by(label_tip_1==sequence_name)) %>%
    rename(x = group) %>%
    left_join(metadata_data,
              by=join_by(label_tip_2==sequence_name)) %>%
    rename(y = group)
  
}

#Make a test duckdb file to tunnel into calculate_rr_matrix
con <- dbConnect(duckdb(),"db_files/db_test_pairs.duckdb",read_only=FALSE)
dbWriteTable(con, "df_pairs", test_pairs, overwrite = TRUE)
df_p <- tbl(con,"df_pairs")

rr <- calculate_rr_matrix(df_p) %>% collect()
dbDisconnect(con)
