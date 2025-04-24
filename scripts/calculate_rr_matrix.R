#File: calculate_rr_matrix.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Function to calculate the identical pair RR matrix, note this version is not optimized and runs in R
#Parameters
#df_p: DB table of pairs of identical sequences with strain names and exposures from bind_pairs_exp
#exp_var: Name of the exposure variable that you want to make the matrix over
#library(doParallel)

calculate_rr_matrix <- function(df_p){

#Assumes the pair matrix is already symmetic
 # df_rr <- df_p %>%
 #   group_by(x,y) %>%
 #   summarise(n_pairs = as.numeric(n())) %>%
 #   group_by(x) %>%
 #   mutate(n_pairs_x = sum(n_pairs)) %>%
 #   group_by(y) %>%
 #   mutate(n_pairs_y = sum(n_pairs)) %>%
 #   ungroup() %>%
 #   mutate(n_pairs_total = sum(n_pairs),
 #          RR = (n_pairs * n_pairs_total + 1) / n_pairs_x / n_pairs_y) %>%
 #   arrange(x,y)
  
  pair_counts <- df_p %>% #Standardize order and count
      group_by(x, y) %>%
      summarise(pair_count = n(), .groups = "drop") %>%
      arrange(x) %>%
      collect()
  
  var_counts <- df_p %>% #Count number of pairs containing each var
    group_by(x) %>%
    summarize(count = n(), .groups = "keep") %>%
    rename(var = x) %>%
    select(var, count) %>%
    collect()

  total_pairs <- df_p %>%
    summarise(total_pairs = n()) %>%
    collect()
  total_pairs <- total_pairs$total_pairs # Convert to a simple numeric

  vars <- var_counts$var

  # Generate all possible symmetric pairs (including self-pairs)
  # Also attaches "join_keys" to make symmetric joining possible
  df_rr <- expand_grid(x = vars, y = vars) %>%
    mutate(
      join_key1 = pmin(x, y), #Standardize how we are joining the counts to the ratio matrix
      join_key2 = pmax(x, y)
    )

  # Join everything and compute relative ratio
  df_rr <- df_rr %>%
    left_join(pair_counts, by = c(join_by(join_key1==x,join_key2==y))) %>%
    mutate(pair_count = replace_na(pair_count, 0)) %>%
    left_join(var_counts, by = c("x" = "var")) %>%
    rename(x_appearances = count) %>%
    left_join(var_counts, by = c("y" = "var")) %>%
    rename(y_appearances = count) %>%
    mutate(
      N_total = total_pairs,
      RR = ((pair_count * N_total) + 1) / (x_appearances * y_appearances)
    ) %>%
    select(x, y, RR, pair_count, x_appearances, y_appearances, N_total) %>%
    rename(x = x) %>%
    rename(y = y) %>%
    rename(N_pairs = pair_count) %>%
    rename(N_x = x_appearances) %>%
    rename(N_y = y_appearances) %>%
    arrange(x)

  return(df_rr)
}
