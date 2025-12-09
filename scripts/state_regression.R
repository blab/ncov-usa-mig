# File: state_regression.R
# Author(s): Amin Bemanian
# Date: 10/2/24
# Description: Comparison of state RRs with interstate travel RRs and regression formation
# Arguments:
#--scenario: Scenario corresponding to data files, for this file only use 50 states + DC geography
# Output: Statistical analysis and regression models comparing different travel metrics with sequence-based RR

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(splines)
  library(broom)
  library(sandwich)
  library(lmtest)
  library(ggplot2)
  library(dplyr)
  library(rsample)
})

# Command line arguments
collect_args <- function() {
  parser <- argparse::ArgumentParser()
  parser$add_argument("--scenario",
    type = "character", default = "CAM_1000",
    help = "Which scenario to perform the analysis on"
  )
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# Read sequence RR data
fn_rr <- paste0("results/", scenario, "/df_RR_by_state.tsv")
state_rr <- fread(fn_rr) %>%
  select(x, y, euclid_dist, min_cbsa_dist, RR) %>%
  rename(RR_seq = RR)

# Load preprocessed travel data
df_travel <- fread("data/travel_vars.tsv") %>%
  inner_join(state_rr, by = c("x", "y")) %>%
  mutate(same_state = x == y) # Add same_state indicator #Join with travel so we have all the values

# Comparison of NHTS and Safegraph RRs
ggplot(df_travel) +
  aes(x = RR_trips, y = RR_move) +
  geom_point(alpha = 0.4, aes(color = same_state)) +
  geom_abline(slope = 1, color = "red", linetype = "dashed") +
  scale_x_continuous(transform = "log10", name = expression(RR["NHTS"]), limits = c(1E-3, 1E3)) +
  scale_y_continuous(transform = "log10", name = expression(RR["Safegraph"]), limits = c(1E-3, 1E3)) +
  theme_bw()
  
ggsave("figs/travel_RR_comparison.jpg", width = 6, height = 5, dpi = 192)

num_states <- nrow(df_travel) %>% sqrt()
k <- 1:(num_states^2)
i <- (k - 1) %% (num_states) + 1
j <- floor(k / num_states) + 1
k_sym <- num_states * (i - 1) + j

# Calculate correlations between sequence RR and mobility metrics
cor_tests <- list(
  NHTS = cor.test(df_travel$RR_seq, df_travel$RR_trips, method = "spearman"),
  Air = cor.test(df_travel$RR_seq, df_travel$RR_air, method = "spearman"),
  Mobility = cor.test(df_travel$RR_seq, df_travel$RR_move,
    method = "spearman",
    use = "complete.obs"
  )
)

# Extract correlation coefficients and p-values
correlations <- data.frame(
  metric = names(cor_tests),
  rho = sapply(cor_tests, function(x) x$estimate),
  p.value = sapply(cor_tests, function(x) x$p.value)
)

message("\nCorrelations with sequence-based RR:")
print(correlations)

# Create output directory
dir.create(file.path("figs", scenario), showWarnings = FALSE)

# Find minimum non-zero values for zero correction
zero_corrections <- list(
  trips = min(df_travel$RR_trips[df_travel$RR_trips > 0], na.rm = TRUE),
  air = min(df_travel$RR_air[df_travel$RR_air > 0], na.rm = TRUE),
  move = min(df_travel$RR_move[df_travel$RR_move > 0], na.rm = TRUE)
)

# Normalize travel RR variables (z-score normalization)
# Distance variables (min_cbsa_dist, euclid_dist) are left unnormalized
normalize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Add normalized travel RR variables to the main df_travel dataframe
df_travel <- df_travel %>%
  mutate(
    zlog_RR_trips = normalize(log10(RR_trips + zero_corrections$trips)),
    zlog_RR_air = normalize(log10(RR_air + zero_corrections$air)),
    zlog_RR_move = normalize(log10(RR_move + zero_corrections$move))
  )

# Add small constant to zeros before taking logs and normalize travel RRs
df_model <- df_travel %>%
  mutate(
    log_RR_seq = log10(RR_seq)
  ) %>%
  filter(!is.infinite(log_RR_seq))

# Fit regression models
models <- list(
  NHTS = lm(log_RR_seq ~ zlog_RR_trips, data = df_model),
  Air = lm(log_RR_seq ~ zlog_RR_air, data = df_model),
  Mobility = lm(log_RR_seq ~ zlog_RR_move, data = df_model),
  Distance = lm(log_RR_seq ~ min_cbsa_dist, data = df_model),
  Combined = lm(log_RR_seq ~ zlog_RR_trips + zlog_RR_air +
    zlog_RR_move + min_cbsa_dist, data = df_model)
)

# Summarize model results
model_summaries <- lapply(models, function(m) {
  coef_data <- tidy(m)
  fit_data <- glance(m)
  list(
    coefficients = coef_data,
    r.squared = fit_data$r.squared,
    adj.r.squared = fit_data$adj.r.squared
  )
})

# Output results
message("\nModel summaries:")
for (name in names(models)) {
  message("\n", name, " model:")
  print(model_summaries[[name]]$coefficients)
  message("Adjusted R-squared: ", round(model_summaries[[name]]$adj.r.squared, 3))
}
OPTIMIZE_AIR <- FALSE
if (OPTIMIZE_AIR) {
  ## Cross validation work by ChatGPT to determine thresholding for air
  cv_score <- function(df, short_cut, long_cut,
                       three_level = TRUE, # FALSE = 2-level
                       K = 10, # k-fold
                       metric = rmse) { # any yardstick fn
    # 1. slice the data -------------------------------------------------
    if (three_level) {
      df_test <- df %>%
        mutate(
          flight_length = case_when(
            min_cbsa_dist < short_cut ~ "Short",
            min_cbsa_dist > long_cut ~ "Long",
            TRUE ~ "Medium"
          )
        ) %>%
        mutate(flight_length = factor(flight_length, levels = c("Short", "Medium", "Long")))
    } else {
      df_test <- df %>%
        mutate(flight_length = case_when(
          min_cbsa_dist < short_cut ~ "Short",
          TRUE ~ "Long"
        )) %>%
        mutate(flight_length = factor(flight_length, levels = c("Short", "Long")))
    }

    # 2. build folds ----------------------------------------------------
    folds <- vfold_cv(df_test, v = K, strata = flight_length)

    # 3. fit & score each fold -----------------------------------------
    fold_stats <- map_dfr(folds$splits, function(split) {
      train <- analysis(split)
      test <- assessment(split)
      mod <- lm(log10(RR_seq) ~ log10(RR_air):flight_length, data = train)
      pred <- tibble(
        truth    = log10(test$RR_seq), # now a real column
        estimate = predict(mod, test)
      )

      tibble(rmse = rmse(pred, truth, estimate)$.estimate)
    })

    mean_rmse <- mean(fold_stats$rmse)
    print(paste0("Short cut: ", short_cut))
    if (three_level) {
      print(paste0("Long cut: ", long_cut))
    }
    print(paste0("Mean RMSE: ", round(mean_rmse, 3)))
    return(mean_rmse)
  }

  # --------------------------------------------------------------------
  # Build the distance grid to search ----------------------------------
  grid <- expand.grid(
    short_cut = seq(100, 2000, by = 25),
    long_cut  = seq(200, 5000, by = 25)
  ) %>%
    filter(short_cut < long_cut)

  # --------------------------------------------------------------------
  # 1) Three-level optimisation ----------------------------------------
  three_results <- grid %>%
    mutate(rmse = pmap_dbl(
      list(short_cut, long_cut),
      cv_score,
      df           = df_travel_air %>% filter(x != y),
      three_level  = TRUE,
      K            = 10
    ))
  ggplot(three_results, aes(x = short_cut, y = long_cut, fill = rmse)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(
      x = "Short Threshold",
      y = "Long Threshold"
    ) +
    theme_bw()
  best3 <- slice_min(three_results, rmse, n = 1)
  best3

  # --------------------------------------------------------------------
  # 2) Two-level optimisation ------------------------------------------
  grid <- expand.grid(
    short_cut = seq(50, 5000, by = 10),
    long_cut = 5000
  )
  two_results <- grid %>%
    mutate(rmse = pmap_dbl(
      list(short_cut, long_cut),
      cv_score,
      df = df_travel_air %>% filter(x != y),
      three_level = FALSE,
      K = 10
    ))
  ggplot(two_results, aes(x = short_cut, y = rmse)) +
    geom_line(aes(color = rmse)) +
    geom_point(aes(color = rmse)) +
    scale_color_viridis_c() +
    labs(x = "Distance Threshold") +
    theme_bw()

  best2 <- slice_min(two_results, rmse, n = 1)
  best2
}

# Set thresholds based on optimization, hardcoding so it does not have to be run each time
df_travel_dist3 <- df_travel %>%
  filter(x != y) %>%
  mutate(flight_length = ifelse(min_cbsa_dist < 300, "Short",
    ifelse(min_cbsa_dist > 1800, "Long", "Medium")
  )) %>%
  mutate(flight_length = factor(flight_length, levels = c("Short", "Medium", "Long")))

# Single threshold version
df_travel_dist2 <- df_travel %>%
  filter(x != y) %>%
  mutate(flight_length = ifelse(min_cbsa_dist < 300, "Short",
    "Long"
  )) %>%
  mutate(flight_length = factor(flight_length, levels = c("Short", "Long")))

df_travel_dist <- df_travel_dist3 # Use three-level for main analysis

travel_comps <- list(
  c("Short", "Medium"),
  c("Medium", "Long")
)

# Plots of RR by travel length
ggplot(df_travel_dist, aes(x = flight_length, y = RR_seq)) +
  geom_boxplot() +
  ggpubr::stat_compare_means( # <-- adds the bars
    comparisons = travel_comps, # list of 2-element vectors
    label = "p.signif", # show “*, **, ns”, etc.
    method = "wilcox.test" # or "t.test", "kruskal.test", …
  ) +
  theme_bw() +
  labs(x = "Travel Length") +
  scale_y_continuous(
    transform = "log",
    name = expression(RR["identical sequences"]),
    breaks = c(0.01, 0.1, 0.5, 1, 2, 10, 100),
    expand = expansion(mult = c(0.18, 0.13)),
    labels = label_number(accuracy = 0.01),
    limits = c(10^(-1.5), 10^(1.5))
  ) +
  geom_hline(yintercept = 1)


ggplot(df_travel_dist, aes(x = flight_length, y = RR_trips)) +
  geom_boxplot() +
  ggpubr::stat_compare_means( # <-- adds the bars
    comparisons = travel_comps, # list of 2-element vectors
    label = "p.signif", # show “*, **, ns”, etc.
    method = "wilcox.test" # or "t.test", "kruskal.test", …
  ) +
  theme_bw() +
  labs(x = "Travel Length") +
  scale_y_continuous(
    transform = "log",
    name = expression(RR["NHTS"]),
    breaks = c(0.01, 0.1, 0.5, 1, 2, 10, 100),
    expand = expansion(mult = c(0.18, 0.13)),
    labels = label_number(accuracy = 0.01),
    limits = c(10^(-1.5), 10^(1.5))
  ) +
  geom_hline(yintercept = 1)

epsilon <- 1E-4

# Make a custom df with rho by flight length
df_air_plot <- df_travel_dist %>%
  group_by(flight_length) %>%
  mutate(
    rho = cor(RR_air + epsilon, RR_seq + epsilon,
      method = "spearman", use = "complete.obs"
    ),
    facet_label = paste0(flight_length, " (rho = ", round(rho, 2), ")")
  ) %>%
  ungroup() %>%
  mutate(facet_label = factor(
    facet_label,
    levels = c(
      paste0("Short (rho = ", round(unique(rho[flight_length == "Short"]), 2), ")"),
      paste0("Medium (rho = ", round(unique(rho[flight_length == "Medium"]), 2), ")"),
      paste0("Long (rho = ", round(unique(rho[flight_length == "Long"]), 2), ")")
    )
  ))

plot_air_length <- ggplot(df_air_plot) +
  geom_point(aes(x = RR_air + epsilon, y = RR_seq + epsilon),
    color = "black", alpha = 0.2
  ) +
  geom_smooth(aes(x = RR_air + epsilon, y = RR_seq + epsilon),
    color = "firebrick", method = "gam", alpha = 0.8
  ) +
  scale_x_continuous(
    transform = "log",
    name = expression(RR["air travel"]),
    breaks = 10^(-2:2),
    labels = c(expression(10^{-2}), expression(10^{-1}), expression(10^{0}), expression(10^{1}), expression(10^{2})), #nolint
    expand = expansion(mult = c(0.18, 0.13)),
    limits = c(10^-2, 10^1)
  ) +
  scale_y_continuous(
    transform = "log",
    name = expression(RR["identical sequences"]),
    breaks = 10^(-2:2),
    labels = c(expression(10^{-2}), #nolint
               expression(10^{-1}), #nolint
               expression(10^{0}), #nolint
               expression(10^{1}), #nolint
               expression(10^{2})),#nolint
    expand = expansion(mult = c(0.18, 0.13)),
    limits = c(10^-1.2, 10^1.2)
  ) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  facet_wrap(vars(facet_label), ncol = 3)
ggsave(plot_air_length,
  filename = paste0("figs/", scenario, "/travel_seq_air_length.jpg"),
  width = 15,
  height = 6,
  dpi = 192
)

zero_correction <- min(df_travel_dist[df_travel_dist$RR_trips > 0, "RR_trips"])

# Model summary function with diagnostics
model_summary <- function(m, calc_bootstrap = TRUE, model_name = NULL, output_file = NULL) {
  # Open output file connection if provided
  if (!is.null(output_file)) {
    sink(output_file, append = TRUE)
    if (!is.null(model_name)) {
      cat(paste0("\n", strrep("=", 80), "\n"))
      cat(paste0("MODEL: ", model_name, "\n"))
      cat(paste0(strrep("=", 80), "\n\n"))
    }
  }

  # Basic model summary
  cat(strrep("-", 80), "\n")
  cat("STANDARD LINEAR MODEL SUMMARY\n")
  cat(strrep("-", 80), "\n")
  print(summary(m))

  # Count predictors and check VIF if multiple predictors
  n_terms <- length(attr(terms(m), "term.labels"))

  if (n_terms >= 2) {
    cat("\n", strrep("-", 80), "\n")
    cat("VARIANCE INFLATION FACTORS (VIF)\n")
    cat(strrep("-", 80), "\n")
    print(car::vif(m))
  } else {
    cat("\n", strrep("-", 80), "\n")
    cat("VARIANCE INFLATION FACTORS (VIF)\n")
    cat(strrep("-", 80), "\n")
    cat("VIF not shown: model has fewer than 2 predictors.\n")
  }

  # Heteroskedasticity test
  cat("\n", strrep("-", 80), "\n")
  cat("BREUSCH-PAGAN TEST FOR HETEROSKEDASTICITY\n")
  cat(strrep("-", 80), "\n")
  print(lmtest::bptest(m))

  # Robust coefficient test with CIs
  cat("\n", strrep("-", 80), "\n")
  cat("ROBUST COEFFICIENT TEST (HC1 ROBUST SEs)\n")
  cat(strrep("-", 80), "\n")
  robust_vcov <- sandwich::vcovHC(m, type = "HC1")
  robust_coef <- lmtest::coeftest(m, vcov = robust_vcov)
  print(robust_coef)

  # Calculate and print robust confidence intervals
  cat("\n95% Confidence Intervals (Robust SEs):\n")
  robust_ci <- confint(robust_coef)
  print(robust_ci)

  # Bootstrap if requested
  if (calc_bootstrap) {
    cat("\n", strrep("-", 80), "\n")
    cat("BOOTSTRAPPED CONFIDENCE INTERVALS FOR COEFFICIENTS\n")
    cat(strrep("-", 80), "\n")
    set.seed(17)  # Set seed for reproducibility
    boot_m <- car::Boot(m, R = 200)
    print(summary(boot_m))
    print(confint(boot_m))
  }

  # AIC
  cat("\n", strrep("-", 80), "\n")
  cat("AKAIKE INFORMATION CRITERION (AIC)\n")
  cat(strrep("-", 80), "\n")
  cat(AIC(m), "\n")

  # Close sink if opened
  if (!is.null(output_file)) {
    sink()
  }
}

# Function to create observed vs predicted plots for RR
# Args:
#   m: Model object
#   RR_obs: Vector of observed RR values
#   model_name: Title for plot
#   log_predict: Whether predictions are log-transformed
#   eps: Small constant for log-safety
plot_pred_RR <- function(m, RR_obs, model_name, log_predict = TRUE, eps = 1e-10) {
  # Generate predictions
  yhat <- predict(m)
  if (log_predict) {
    yhat <- 10^yhat # Convert log predictions back to original scale
  }

  # Prepare data for plotting
  df <- data.frame(
    RR_obs = as.numeric(RR_obs),
    RR_hat = as.numeric(yhat)
  )
  df <- df[complete.cases(df), ] # Remove NA values
  df$RR_obs <- df$RR_obs + eps
  df$RR_hat <- df$RR_hat + eps

  # Calculate RMSE
  rmse_val <- yardstick::rmse_vec(truth = df$RR_obs, estimate = df$RR_hat)
  rmse_val <- as.numeric(rmse_val)
  message(sprintf("RMSE: %.3f", rmse_val))

  # Set plot range
  range_vals <- range(c(df$RR_obs, df$RR_hat), na.rm = TRUE)

  # Create plot
  p <- ggplot(data = df) +
    aes(x = .data$RR_hat, y = .data$RR_obs) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(transform = "log10", name = "Predicted") +
    scale_y_continuous(transform = "log10", name = "Observed") +
    ggtitle(sprintf("%s - RMSE: %.3f", model_name, rmse_val)) +
    coord_equal(xlim = range_vals, ylim = range_vals) +
    theme_bw()

  return(p)
}


df_model_data <- df_travel_dist %>%
  filter(x != y) %>%
  na.omit()  # Remove all rows with any NA values for car::Boot compatibility
y_obs <- df_model_data$RR_seq

# Create output file for regression results
results_file <- paste0("results/", scenario, "/regression_results.txt")
# Clear file if it exists
if (file.exists(results_file)) file.remove(results_file)

lm_trips <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_trips + zero_corrections$trips)
)
lm_trips %>% model_summary(model_name = "NHTS Trips", output_file = results_file)
lm_trips_plot <- plot_pred_RR(lm_trips, y_obs, "NHTS Trips")

slm_trips <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ bs(log10(RR_trips + zero_corrections$trips), df = 6)
)
slm_trips %>% model_summary(model_name = "NHTS Trips - Spline", output_file = results_file)
slm_trips_plot <- plot_pred_RR(slm_trips, y_obs, "NHTS Trips - Spline")

lm_move <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_move + zero_corrections$move)
)
lm_move %>% model_summary(model_name = "Safegraph Mobility", output_file = results_file)
lm_move_plot <- plot_pred_RR(lm_move, y_obs, "Safegraph Mobility")

slm_move <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ bs(log10(RR_move + zero_corrections$move), df = 6)
)
slm_move %>% model_summary(calc_bootstrap = FALSE, model_name = "Safegraph Mobility - Spline", output_file = results_file)
slm_move_plot <- plot_pred_RR(slm_move, y_obs, "Safegraph Mobility - Spline")

lm_air <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_air + zero_corrections$air)
)
lm_air %>% model_summary(model_name = "Air Travel", output_file = results_file)
lm_air_plot <- plot_pred_RR(lm_air, y_obs, "Air Travel")

slm_air <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ bs(log10(RR_air + zero_corrections$air), df = 6)
)
slm_air %>% model_summary(calc_bootstrap = FALSE, model_name = "Air Travel - Spline", output_file = results_file)
slm_air_plot <- plot_pred_RR(slm_air, y_obs, "Air Travel - Spline")

lm_air_int <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_air + zero_corrections$air):flight_length
)
lm_air_int %>% model_summary(model_name = "Air Travel : Flight Length", output_file = results_file)
lm_air_int_plot <- plot_pred_RR(lm_air_int, y_obs, "Air Travel : Flight Length")


slm_air_int <- lm(
  data = df_model_data,
  log10(RR_seq) ~ bs(log10(RR_air + zero_corrections$air), df = 6) * bs(min_cbsa_dist, df = 6)
)
slm_air_int %>% model_summary(calc_bootstrap = FALSE, model_name = "Air Travel x CBSA Distance - Spline", output_file = results_file)
slm_air_int_plot <- plot_pred_RR(slm_air_int, y_obs, "Air Travel x CBSA Distance - Spline")

lm_euclid <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ euclid_dist
)
lm_euclid %>% model_summary(model_name = "Euclidean Distance", output_file = results_file)
lm_euclid_plot <- plot_pred_RR(lm_euclid, y_obs, "Euclidean Distance")

slm_euclid <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ bs(euclid_dist, df = 6)
)
slm_euclid %>% model_summary(calc_bootstrap = FALSE, model_name = "Euclidean Distance - Spline", output_file = results_file)
slm_euclid_plot <- plot_pred_RR(slm_euclid, y_obs, "Euclidean Distance - Spline")

lm_cbsa <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ min_cbsa_dist
)
lm_cbsa %>% model_summary(model_name = "CBSA Distance", output_file = results_file)
lm_cbsa_plot <- plot_pred_RR(lm_cbsa, y_obs, "CBSA Distance")

slm_cbsa <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ bs(min_cbsa_dist, df = 6)
)
slm_cbsa %>% model_summary(calc_bootstrap = FALSE, model_name = "CBSA Distance - Spline", output_file = results_file)
slm_cbsa_plot <- plot_pred_RR(slm_cbsa, y_obs, "CBSA Distance - Spline")

# Combined Regression
lm_combo <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ zlog_RR_trips +
    zlog_RR_move +
    zlog_RR_air:flight_length
)
lm_combo %>% model_summary(model_name = "Multivariate Model (normalized)", output_file = results_file)
lm_combo_plot <- plot_pred_RR(lm_combo, y_obs, "Multivariate Model (normalized)")

# ============================================================================
# ANOMALY DETECTION FOR STATE PAIRS
# ============================================================================
message("\n", strrep("=", 80))
message("ANOMALY DETECTION: State pairs enriched beyond model expectations")
message(strrep("=", 80), "\n")

# Calculate diagnostic statistics
df_anomalies <- df_model_data %>%
  mutate(
    fitted = predict(lm_combo),
    residual = log10(RR_seq) - fitted,
    standardized_residual = rstandard(lm_combo),
    studentized_residual = rstudent(lm_combo),
    cooks_distance = cooks.distance(lm_combo),
    leverage = hatvalues(lm_combo),
    predicted_RR = 10^fitted
  )

# Define anomaly thresholds
STUDENTIZED_THRESHOLD <- 3  # |studentized residual| > 3
COOKS_THRESHOLD <- 4 / nrow(df_model_data)  # Cook's D > 4/n (common rule)

# Flag anomalies
df_anomalies <- df_anomalies %>%
  mutate(
    is_outlier_residual = abs(studentized_residual) > STUDENTIZED_THRESHOLD,
    is_influential = cooks_distance > COOKS_THRESHOLD,
    is_anomaly = is_outlier_residual | is_influential,
    anomaly_type = case_when(
      is_outlier_residual & residual > 0 ~ "Enriched",
      is_outlier_residual & residual < 0 ~ "Depleted",
      is_influential ~ "Influential",
      TRUE ~ "Normal"
    )
  ) %>%
  arrange(desc(abs(studentized_residual)))

# Summary statistics
n_enriched <- sum(df_anomalies$anomaly_type == "Enriched")
n_depleted <- sum(df_anomalies$anomaly_type == "Depleted")
n_influential <- sum(df_anomalies$is_influential & !df_anomalies$is_outlier_residual)

message(sprintf("Total state pairs analyzed: %d", nrow(df_anomalies)))
message(sprintf("Enriched pairs (higher than expected): %d (%.1f%%)",
                n_enriched, 100 * n_enriched / nrow(df_anomalies)))
message(sprintf("Depleted pairs (lower than expected): %d (%.1f%%)",
                n_depleted, 100 * n_depleted / nrow(df_anomalies)))
message(sprintf("Influential pairs: %d (%.1f%%)",
                n_influential, 100 * n_influential / nrow(df_anomalies)))

# Top enriched pairs (transmission higher than expected from travel)
df_enriched <- df_anomalies %>%
  filter(anomaly_type == "Enriched") %>%
  select(x, y, RR_seq, predicted_RR, residual, studentized_residual,
         cooks_distance, RR_trips, RR_air, RR_move, min_cbsa_dist) %>%
  arrange(desc(studentized_residual))

message("\nTop enriched state pairs (observed > expected):")
print(head(df_enriched, 20))

# Top depleted pairs (transmission lower than expected from travel)
df_depleted <- df_anomalies %>%
  filter(anomaly_type == "Depleted") %>%
  select(x, y, RR_seq, predicted_RR, residual, studentized_residual,
         cooks_distance, RR_trips, RR_air, RR_move, min_cbsa_dist) %>%
  arrange(studentized_residual)

message("\nTop depleted state pairs (observed < expected):")
print(head(df_depleted, 20))

# Save anomaly results
anomaly_file <- paste0("results/", scenario, "/state_pair_anomalies.tsv")
fwrite(df_anomalies, anomaly_file, sep = "\t")
message(sprintf("\nAnomaly results saved to: %s", anomaly_file))

# ============================================================================
# VISUALIZATIONS
# ============================================================================

# 1. Residual plot with anomalies highlighted
# Remove symmetric pairs (keep only alphabetically ordered pairs)
df_anomalies_unique <- df_anomalies %>%
  mutate(
    edge_pair = paste(pmin(x, y), pmax(x, y), sep = "-"),
    pair_label = paste0(pmin(x, y), "-", pmax(x, y))
  ) %>%
  group_by(edge_pair) %>%
  slice_head(n = 1) %>%  # Keep first occurrence of each unique pair
  ungroup() %>%
  mutate(
    label = ifelse(anomaly_type %in% c("Enriched", "Depleted"), pair_label, "")
  )

# Update counts for subtitle after deduplication
n_enriched_unique <- sum(df_anomalies_unique$anomaly_type == "Enriched")
n_depleted_unique <- sum(df_anomalies_unique$anomaly_type == "Depleted")

plot_residuals <- ggplot(df_anomalies_unique, aes(x = fitted, y = residual)) +
  geom_point(aes(color = anomaly_type, size = cooks_distance), alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = c(-STUDENTIZED_THRESHOLD, STUDENTIZED_THRESHOLD) * sd(df_anomalies_unique$residual),
             linetype = "dashed", color = "red", alpha = 0.5) +
  ggrepel::geom_text_repel(
    aes(label = label, color = anomaly_type),
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.2,
    segment.alpha = 0.6,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Enriched" = "firebrick", "Depleted" = "steelblue",
               "Influential" = "orange", "Normal" = "gray60"),
    breaks = c("Enriched", "Depleted", "Influential", "Normal")
  ) +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(
    x = "Fitted log10(RR)",
    y = "Residual",
    title = "State Pair Anomalies in Multivariate Model",
    subtitle = paste0("Enriched: n=", n_enriched_unique, ", Depleted: n=", n_depleted_unique, " (unique pairs)"),
    color = "Anomaly Type",
    size = "Cook's Distance"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave(paste0("figs/", scenario, "/state_pair_anomalies_residuals.jpg"),
       plot = plot_residuals, width = 10, height = 7, dpi = 192)

# 2. Cook's Distance plot
plot_cooks <- ggplot(df_anomalies, aes(x = 1:nrow(df_anomalies), y = cooks_distance)) +
  geom_segment(aes(xend = 1:nrow(df_anomalies), yend = 0, color = is_influential),
               alpha = 0.6) +
  geom_point(aes(color = is_influential), size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = COOKS_THRESHOLD, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick")) +
  labs(
    x = "Observation Index",
    y = "Cook's Distance",
    title = "Influential State Pairs (Cook's Distance)",
    subtitle = sprintf("Threshold: 4/n = %.4f", COOKS_THRESHOLD),
    color = "Influential"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(paste0("figs/", scenario, "/state_pair_cooks_distance.jpg"),
       plot = plot_cooks, width = 10, height = 6, dpi = 192)

# 3. Q-Q plot for residuals
plot_qq <- ggplot(df_anomalies, aes(sample = studentized_residual)) +
  stat_qq(aes(color = is_anomaly), alpha = 0.6) +
  stat_qq_line(color = "black", linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick")) +
  labs(
    title = "Q-Q Plot: Studentized Residuals",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles",
    color = "Anomaly"
  ) +
  theme_bw()

ggsave(paste0("figs/", scenario, "/state_pair_qq_plot.jpg"),
       plot = plot_qq, width = 8, height = 6, dpi = 192)

# 4. Observed vs Predicted with anomalies highlighted
plot_obs_pred_anomaly <- ggplot(df_anomalies, aes(x = predicted_RR, y = RR_seq)) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_point(aes(color = anomaly_type, size = abs(studentized_residual)), alpha = 0.6) +
  scale_x_continuous(transform = "log10", name = "Predicted RR") +
  scale_y_continuous(transform = "log10", name = "Observed RR") +
  scale_color_manual(
    values = c("Enriched" = "firebrick", "Depleted" = "steelblue",
               "Influential" = "orange", "Normal" = "gray60")
  ) +
  scale_size_continuous(range = c(0.5, 4), name = "|Studentized Residual|") +
  coord_equal() +
  labs(
    title = "Observed vs Predicted RR with Anomalies",
    color = "Anomaly Type"
  ) +
  theme_bw()

ggsave(paste0("figs/", scenario, "/state_pair_obs_pred_anomalies.jpg"),
       plot = plot_obs_pred_anomaly, width = 10, height = 8, dpi = 192)

message("\nAnomaly visualizations saved to figs/", scenario, "/")

# ============================================================================
# ROBUST REGRESSION: REFIT MODEL WITHOUT ANOMALIES
# ============================================================================
message("\n", strrep("=", 80))
message("ROBUST MODEL: Refitting lm_combo without anomalies/influential obs")
message(strrep("=", 80), "\n")

# Filter out anomalies for robust model
df_model_robust <- df_model_data %>%
  filter(!df_anomalies$is_anomaly)

y_obs_robust <- df_model_robust$RR_seq

n_removed <- nrow(df_model_data) - nrow(df_model_robust)
pct_removed <- 100 * n_removed / nrow(df_model_data)

message(sprintf("Original sample size: %d", nrow(df_model_data)))
message(sprintf("Removed observations: %d (%.1f%%)", n_removed, pct_removed))
message(sprintf("Robust sample size: %d\n", nrow(df_model_robust)))

# Refit the combined model without anomalies
lm_combo_robust <- lm(
  data = df_model_robust,
  formula = log10(RR_seq) ~ zlog_RR_trips +
    zlog_RR_move +
    zlog_RR_air:flight_length
)

lm_combo_robust %>% model_summary(model_name = "Multivariate Model - Robust (anomalies removed)",
                                   output_file = results_file)
lm_combo_robust_plot <- plot_pred_RR(lm_combo_robust, y_obs_robust,
                                       "Multivariate Model - Robust (anomalies removed)")

# ============================================================================
# MODEL COMPARISON: ORIGINAL vs ROBUST
# ============================================================================
message("\n", strrep("=", 80))
message("MODEL COMPARISON: Original vs Robust")
message(strrep("=", 80), "\n")

# Extract key statistics for comparison
comparison_stats <- data.frame(
  Model = c("Original (with anomalies)", "Robust (anomalies removed)"),
  N = c(nrow(df_model_data), nrow(df_model_robust)),
  R_squared = c(summary(lm_combo)$r.squared, summary(lm_combo_robust)$r.squared),
  Adj_R_squared = c(summary(lm_combo)$adj.r.squared, summary(lm_combo_robust)$adj.r.squared),
  RMSE = c(
    sqrt(mean(residuals(lm_combo)^2)),
    sqrt(mean(residuals(lm_combo_robust)^2))
  ),
  AIC = c(AIC(lm_combo), AIC(lm_combo_robust))
)

print(comparison_stats)

# Coefficient comparison
coef_comparison <- data.frame(
  Term = names(coef(lm_combo)),
  Original = coef(lm_combo),
  Robust = c(coef(lm_combo_robust), rep(NA, length(coef(lm_combo)) - length(coef(lm_combo_robust))))
)
# Handle case where models have same terms
if (length(coef(lm_combo)) == length(coef(lm_combo_robust))) {
  coef_comparison$Robust <- coef(lm_combo_robust)
  coef_comparison$Difference <- coef_comparison$Robust - coef_comparison$Original
  coef_comparison$Pct_Change <- 100 * coef_comparison$Difference / abs(coef_comparison$Original)
}

message("\nCoefficient comparison:")
print(coef_comparison)

# Save comparison stats
comparison_file <- paste0("results/", scenario, "/model_comparison_original_vs_robust.tsv")
fwrite(comparison_stats, comparison_file, sep = "\t")
message(sprintf("\nModel comparison saved to: %s", comparison_file))

# ============================================================================
# VISUALIZATION: SIDE-BY-SIDE COMPARISON
# ============================================================================

# Create combined plot showing both models
library(patchwork)
comparison_plot <- lm_combo_plot + lm_combo_robust_plot +
  plot_annotation(
    title = "Model Comparison: With vs Without Anomalies/Influential Obs",
    subtitle = sprintf("Removed %d anomalies/influential obs (%.1f%% of data)", n_removed, pct_removed)
  )

ggsave(paste0("figs/", scenario, "/model_comparison_with_without_anomalies.jpg"),
       plot = comparison_plot, width = 14, height = 6, dpi = 192)

message("\nComparison plot saved to figs/", scenario, "/model_comparison_with_without_anomalies.jpg")

# Spline combined regression
slm_combo <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ zlog_RR_trips +
    zlog_RR_move +
    bs(zlog_RR_air, df = 6) * bs(min_cbsa_dist, df = 6)
)
slm_combo %>% model_summary(calc_bootstrap = FALSE, model_name = "Multivariate Model - Spline (normalized)", output_file = results_file)
slm_combo_plot <- plot_pred_RR(slm_combo, y_obs, "Multivariate Model - Spline (normalized)")

OE_plots <- patchwork::wrap_plots(lm_trips_plot, slm_trips_plot, lm_move_plot, slm_move_plot,
  lm_air_plot, slm_air_plot, lm_air_int_plot, slm_air_int_plot,
  lm_euclid_plot, slm_euclid_plot, lm_cbsa_plot, slm_cbsa_plot,
  patchwork::plot_spacer(), lm_combo_plot, slm_combo_plot, patchwork::plot_spacer(),
  ncol = 4, nrow = 4
)
ggsave(paste0("figs/", scenario, "/dist/regression_fit.jpg"),
  plot = OE_plots,
  width = 22,
  height = 20,
  units = "in",
  dpi = 300
)
