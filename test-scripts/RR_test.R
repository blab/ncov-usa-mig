#File: RR_test.R
#Author(s): Amin Bemanian
#Date: 4/20/24
#Description: Creates a test database for debugging the RR calculaltion code

library(tidyverse)
library(data.table)
library(dbplyr)
library(duckdb)

source("scripts/calculate_rr_matrix.R")

#Make the test table
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

#Make a test duckdb file to tunnel into calculate_rr_matrix
con <- dbConnect(duckdb(),"db_files/db_test_pairs.duckdb",read_only=FALSE)
dbWriteTable(con, "df_pairs", test_pairs, overwrite = TRUE)
df_p <- tbl(con,"df_pairs")

calculate_rr_matrix(df_p)
dbDisconnect(con)
