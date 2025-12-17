#File: calculate_rr_matrix.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Function to calculate the identical pair RR matrix, uses DuckDB
#Parameters
#df_p: DB table of pairs of identical sequences with strain names and exposures from bind_pairs_exp
#exp_var: Name of the exposure variable that you want to make the matrix over
#library(doParallel)

calculate_rr_matrix <- function(df_p){

  DEBUG_MODE <- FALSE #Use for debugging, by design not provided as an argument

  if(DEBUG_MODE){
    # Check if df_p is already symmetric
    # For each pair (strain_1=A, strain_2=B), check if reverse pair (strain_1=B, strain_2=A) exists

    # Get all distinct off-diagonal pairs
    forward_pairs <- df_p %>%
      filter(strain_1 != strain_2) %>%
      distinct(strain_1, strain_2) %>%
      summarise(n = n()) %>%
      pull(n)

    # Check how many of these have their reverse
    # Join forward pairs (A,B) with all pairs to find (B,A)
    reversed_pairs <- df_p %>%
      filter(strain_1 != strain_2) %>%
      select(s1 = strain_1, s2 = strain_2) %>%
      distinct() %>%
      inner_join(
        df_p %>% select(strain_1, strain_2) %>% distinct(),
        by = c("s2" = "strain_1", "s1" = "strain_2")
      ) %>%
      summarise(n = n()) %>%
      pull(n)

    is_symmetric <- (forward_pairs == reversed_pairs && forward_pairs > 0)
    message(sprintf("Pairs symmetric check: %d forward pairs, %d have reverses. Symmetric: %s",
                    forward_pairs, reversed_pairs, is_symmetric))
  }

  # Default is to assume the pairs provided are not symmetric so we symmetrize it
  # Symmetrize df_p: add reverse of each off-diagonal pair
  # This is needed because pairs table stores each identical pair once (A,B)
  # but RR calculation conceptually needs both directions (A→B and B→A)
  df_p <- df_p %>%
    union_all(
      df_p %>%
        mutate(
          x_temp = x, y_temp = y,
          strain_1_temp = strain_1, strain_2_temp = strain_2
        ) %>%
        mutate(
          x = y_temp, y = x_temp,
          strain_1 = strain_2_temp, strain_2 = strain_1_temp
        ) %>%
        select(-x_temp, -y_temp, -strain_1_temp, -strain_2_temp)
    )

  # Count number of pairs for each x and y combination
  pair_counts <- df_p %>%
      group_by(x, y) %>%
      summarise(pair_count = n(), .groups = "drop") %>%
      arrange(x, y) %>%
      collect()
  if(DEBUG_MODE){
    print("Sum of pair counts:")
    print(sum(pair_counts$pair_count))
  }

  #Count number of pairs containing each var
  var_counts <- df_p %>%
    group_by(x) %>%
    summarize(count = n(), .groups = "keep") %>%
    rename(var = x) %>%
    select(var, count) %>%
    arrange(var) %>%
    collect()

  if(DEBUG_MODE){
    print("List of var values")
    print(var_counts)
    print("Total var counts")
    print(sum(var_counts$count))
  }

  total_pairs <- df_p %>%
    summarise(total_pairs = n()) %>%
    collect()
  total_pairs <- total_pairs$total_pairs # Convert to a simple numeric
  if(DEBUG_MODE){
    print("Total pair length:")
    print(total_pairs)
  }
  vars <- var_counts$var

  # Generate all possible pairs (including self-pairs)
  df_rr <- expand_grid(x = vars, y = vars)

  # Join everything and compute relative ratio
  df_rr <- df_rr %>%
    left_join(pair_counts, by = c("x", "y")) %>%
    mutate(pair_count = replace_na(pair_count, 0)) %>%
    left_join(var_counts, by = c("x" = "var")) %>%
    rename(x_appearances = count) %>%
    mutate(x_appearances = as.numeric(x_appearances)) %>%
    left_join(var_counts, by = c("y" = "var")) %>%
    rename(y_appearances = count) %>%
    mutate(y_appearances = as.numeric(y_appearances)) %>%
    mutate(
      N_total = as.numeric(total_pairs),
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
