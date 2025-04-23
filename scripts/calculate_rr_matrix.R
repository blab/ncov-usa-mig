#File: calculate_rr_matrix.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Function to calculate the identical pair RR matrix, note this version is not optimized and runs in R
#Parameters
#df_p: DB table of pairs of identical sequences with strain names and exposures from bind_pairs_exp
#exp_var: Name of the exposure variable that you want to make the matrix over
#library(doParallel)

calculate_rr_matrix <- function(df_p){

#df_p_symmetric <- df_p %>% #Mirror matrix with the reverse ordering of X and Y
#  rename(x_future = y, y = x) %>%
#  rename(x = x_future)

#  Old version entirely in DuckDB
#  df_rr <- df_p_symmetric %>%
#    union_all(df_p) %>%
#    group_by(x,y) %>%
#    summarise(n_pairs = as.numeric(n())) %>%
#    group_by(x) %>%
#    mutate(n_pairs_x = sum(n_pairs)) %>%
#    group_by(y) %>%
#    mutate(n_pairs_y = sum(n_pairs)) %>%
#    ungroup() %>%
#    mutate(n_pairs_total = sum(n_pairs),
#           RR = (n_pairs * n_pairs_total + 1) / n_pairs_x / n_pairs_y) %>%
#    arrange(x,y)
  
  pair_counts <- df_p %>% #Standardize order and count
      mutate(
        x1 = if_else(x < y, x, y),
        x2 = if_else(x >= y, x, y)
      ) %>%
      group_by(x1, x2) %>%
      summarise(pair_count = n(), .groups = "drop") %>%
      arrange(x1) %>%
      collect()

  var_counts <- df_p %>% #Count number of pairs containing each var
    mutate(rowid = row_number()) %>%
    select(rowid, x, y) %>%
    pivot_longer(cols = c(x, y), names_to = "which", values_to = "var") %>%
    distinct(rowid, var) %>%
    group_by(var) %>%
    summarise(pair_appearances = n(), .groups = "drop") %>%
    arrange(var) %>%
    collect()

  total_pairs <- df_p %>% 
    summarise(total_pairs = n()) %>%
    collect()
  total_pairs <- total_pairs$total_pairs # Convert to a simple numeric
  
  vars <- var_counts$var
  
  # Generate all possible symmetric pairs (including self-pairs)
  # Also attaches "join_keys" to make symmetric joining possible
  df_rr <- expand_grid(x1 = vars, x2 = vars) %>%
    mutate(
      join_key1 = pmin(x1, x2), #Standardize how we are joining the counts to the ratio matrix
      join_key2 = pmax(x1, x2)
    )
  
  # Join everything and compute relative ratio
  df_rr <- df_rr %>%
    left_join(pair_counts, by = c(join_by(join_key1==x1,join_key2==x2))) %>%
    mutate(pair_count = replace_na(pair_count, 0)) %>%
    left_join(var_counts, by = c("x1" = "var")) %>%
    rename(x1_appearances = pair_appearances) %>%
    left_join(var_counts, by = c("x2" = "var")) %>%
    rename(x2_appearances = pair_appearances) %>%
    mutate(
      N_total = total_pairs,
      RR = ((pair_count * N_total) + 1) / (x1_appearances * x2_appearances)
    ) %>%
    select(x1, x2, RR, pair_count, x1_appearances, x2_appearances, N_total) %>%
    rename(x = x1) %>%
    rename(y = x2) %>%
    rename(N_pairs = pair_count) %>%
    rename(N_x = x1_appearances) %>%
    rename(N_y = x2_appearances) %>%
    arrange(x)

  return(df_rr)
}
