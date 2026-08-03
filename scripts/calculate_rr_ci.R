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

source("scripts/subsample_ci.R")

calculate_rr_ci <- function(db_con, exp_var, samp_cov = 0.8, k = 200,
                            interval_width = 0.95, exclude_duplicates = TRUE) {

  rr_list <- vector("list", k)

  for (i in 1:k) {
    tictoc::tic()
    iteration_pairs <- bind_pairs_exp(db_con, exp_var, sub_samp = TRUE, samp_cov,
                                   exclude_duplicates = exclude_duplicates)
    iteration_rr <- iteration_pairs %>% 
      calculate_rr_matrix() %>%
      collect()
    iteration_rr_cleaned <- iteration_rr %>% 
      select(x, y, RR) %>%
      rename(!!paste0("RR_", i) := RR)
    rr_list[[i]] <- iteration_rr_cleaned
    if (i %% 10 == 0) {
      gc()
      message(sprintf("Completed bootstrap iteration %d/%d", i, k))
    }
    tictoc::toc()
  }


  # Join all iterations by (x, y) key
  rr_combined <- reduce(rr_list, full_join, by = c("x", "y"))

  # Rescaled quantiles across columns (see scripts/subsample_ci.R for why the
  # raw replicate spread would be too narrow by sqrt((1 - samp_cov)/samp_cov)).
  rr_matrix <- rr_combined %>% select(starts_with("RR_")) %>% as.matrix()
  ci_bounds <- t(apply(rr_matrix, 1, subsample_ci_bounds,
                       samp_cov = samp_cov, interval_width = interval_width))

  rr_combined %>%
    mutate(
      ci_lb = ci_bounds[, 1],
      ci_ub = ci_bounds[, 2]
    ) %>%
    select(x, y, ci_lb, ci_ub)
}
