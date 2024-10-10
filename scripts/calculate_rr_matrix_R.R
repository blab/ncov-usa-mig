#File: calculate_rr_matrix.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Function to calculate the identical pair RR matrix, note this version is not optimized and runs in R
#Parameters
#df_p: Dataframe of pairs of identical sequences with strain names and exposures from bind_pairs_exp
#exp_var: Name of the exposure variable that you want to make the matrix over
library(doParallel)

calculate_rr_matrix <- function(df_p,exp_var,
                                return_ci = FALSE,
                                prop_subsample = 0.8,
                                n_subsample = 1000,
                                interval_width = 0.95){

  
  df_p$x <- df_p[[paste0("X_",exp_var)]]
  df_p$y <- df_p[[paste0("Y_",exp_var)]]
  df_p_symmetric <- df_p %>% #Mirror matrix with the reverse ordering of X and Y
    rename(x_future = y, y = x) %>%
    rename(x = x_future)
  df_rr <- df_p_symmetric %>% 
    bind_rows(df_p) %>% 
    group_by(x,y) %>% 
    summarise(n_pairs = as.numeric(n())) %>% 
    group_by(x) %>% 
    mutate(n_pairs_x = sum(n_pairs)) %>% 
    group_by(y) %>% 
    mutate(n_pairs_y = sum(n_pairs)) %>% 
    ungroup() %>% 
    mutate(n_pairs_total = sum(n_pairs),
           RR = n_pairs  * n_pairs_total / n_pairs_x / n_pairs_y)
  
  if(return_ci == TRUE){ #Returning a credible interval for RR
    rr_dist <- foreach(i=1:n_subsample,.combine = cbind) %dopar% {
      df_samp <- sample_frac(df_p,prop_subsample) #Subsample
      rr_samp <- calculate_rr_matrix(df_samp,exp_var) #Recursively run on the subsample
      rr_samp$RR #Return the relative risk as a wide distribution
    }
    ci_lb <- apply(rr_dist,1,quantile,probs=(1-interval_width)/2)
    ci_ub <- apply(rr_dist,1,quantile,probs=(1+interval_width)/2)
    df_rr <- cbind(df_rr,ci_lb,ci_ub)
  }
  
  return(df_rr)
}