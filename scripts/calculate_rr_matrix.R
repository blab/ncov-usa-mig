#File: calculate_rr_matrix.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Function to calculate the identical pair RR matrix, note this version is not optimized and runs in R
#Parameters
#df_p: DB table of pairs of identical sequences with strain names and exposures from bind_pairs_exp
#exp_var: Name of the exposure variable that you want to make the matrix over
#library(doParallel)

calculate_rr_matrix <- function(df_p){

  df_p_symmetric <- df_p %>% #Mirror matrix with the reverse ordering of X and Y
    rename(x_future = y, y = x) %>%
    rename(x = x_future)
  df_rr <- df_p_symmetric %>% 
    union_all(df_p) %>% 
    group_by(x,y) %>% 
    summarise(n_pairs = as.numeric(n())) %>% 
    group_by(x) %>% 
    mutate(n_pairs_x = sum(n_pairs)) %>% 
    group_by(y) %>% 
    mutate(n_pairs_y = sum(n_pairs)) %>% 
    ungroup() %>% 
    mutate(n_pairs_total = sum(n_pairs),
           RR = n_pairs  * n_pairs_total / n_pairs_x / n_pairs_y) %>%
    arrange(x,y)
  
  return(df_rr)
}
