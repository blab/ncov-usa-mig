#File: calculate_rr_ci.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Wrapper for RR confidence interval calculation
#works by subsetting relies on bind_pairs_exp
#TO-DO: Add time series support (time_bounds)
#Arguments: None
#db_con: DuckDB Connection
#exp_var: Name of exposure variable to attach
#sub_samp: Sub-sample (for CI calculation)
#samp_cov: Sample coverage (for CI calculation)
#k: Number of repeated sub-samples
#interval_width: CI interval width
#time_bounds: Date range (as vector of two Date values) for time analyses (NEEDS TO BE ADDED)
#exclude_duplicates: If TRUE, exclude pairs marked as possible_duplicates (default: FALSE)

calculate_rr_ci <- function(db_con,exp_var,samp_cov = 0.8, k = 200, interval_width = 0.95, exclude_duplicates = FALSE){
  RR_dist <- matrix()
  for(x in 1:k){
    RR_samp <- bind_pairs_exp(db_con,exp_var,sub_samp = TRUE,samp_cov, exclude_duplicates = exclude_duplicates) %>%
      calculate_rr_matrix %>%
      collect
    if(x == 1){
      CI_keys <- RR_samp %>% select(c(x,y))
      RR_dist <- RR_samp$RR
    }else{
      RR_dist <- cbind(RR_dist,RR_samp$RR)
    }
  }
  ci_lb <- apply(RR_dist,1,quantile,probs=(1-interval_width)/2)
  ci_ub <- apply(RR_dist,1,quantile,probs=(1+interval_width)/2)
  return(cbind(CI_keys,ci_lb,ci_ub))
}
