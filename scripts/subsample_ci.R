#File: subsample_ci.R
#Author(s): Amin Bemanian
#Description: Shared helpers for turning subsampling replicates into confidence
# intervals, with the finite-sample rescaling applied.
#
#Why subsampling rather than a with-replacement bootstrap: a strain drawn twice
#would join against `pairs` and manufacture a spurious identical pair of the
#strain with its own copy, inflating counts exactly on the heatmap diagonal.
#Drawing a fraction f of strains WITHOUT replacement avoids that entirely.
#
#Why the rescaling: for m-out-of-n sampling without replacement,
#  Var(theta_m) ~ (sigma^2 / m) * (1 - m/n)
#while the quantity of interest is Var(theta_n) ~ sigma^2 / n. Raw quantiles of
#the replicates therefore understate the sampling variability by a factor of
#sqrt((1 - f)/f) -- a factor of 2 at the f = 0.8 used throughout this project.
#We correct by expanding each replicate's deviation from the replicate median by
#  c = sqrt(f / (1 - f))
#which is 1 at f = 0.5 (where no correction is needed) and 2 at f = 0.8.
#
#Expansion is done on the log scale: RR is strictly positive and right-skewed,
#and every RR analysis here is read on a log scale, so this keeps the lower
#bound positive and the expansion symmetric in the space actually plotted.

# Multiplier applied to log-scale deviations from the replicate median.
subsample_ci_scale <- function(samp_cov){
  stopifnot(samp_cov > 0, samp_cov < 1)
  sqrt(samp_cov / (1 - samp_cov))
}

# x: replicate estimates for a single cell (strictly positive)
# Returns c(lb, ub); NA if there are too few usable replicates.
subsample_ci_bounds <- function(x, samp_cov = 0.8, interval_width = 0.95){
  x <- x[is.finite(x) & x > 0]
  if(length(x) < 2){
    return(c(NA_real_, NA_real_))
  }
  log_x <- log(x)
  centre <- median(log_x)
  probs <- c((1 - interval_width) / 2, (1 + interval_width) / 2)
  q <- quantile(log_x, probs = probs, na.rm = TRUE)
  unname(exp(centre + subsample_ci_scale(samp_cov) * (q - centre)))
}

# Scalar wrappers, convenient inside dplyr::summarize().
subsample_ci_lb <- function(x, samp_cov = 0.8, interval_width = 0.95){
  subsample_ci_bounds(x, samp_cov, interval_width)[1]
}

subsample_ci_ub <- function(x, samp_cov = 0.8, interval_width = 0.95){
  subsample_ci_bounds(x, samp_cov, interval_width)[2]
}
