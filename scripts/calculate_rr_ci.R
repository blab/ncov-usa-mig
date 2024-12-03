calculate_rr_ci <- function(db_con,exp_var,samp_cov = 0.8, k = 200, interval_width = 0.95){
  RR_dist <- matrix()
  for(x in 1:k){
    RR_samp <- bind_pairs_exp(db_con,exp_var,sub_samp = TRUE,samp_cov) %>%
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
