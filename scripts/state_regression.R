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
  parser$add_argument("--splines",
    type = "logical", default = FALSE,
    help = "Fit spline variants of each model (default: FALSE)"
  )
  return(parser$parse_args())
}

args <- collect_args()
scenario    <- args$scenario
spline_flag <- args$splines

# Read sequence RR data
fn_rr <- paste0("results/", scenario, "/df_RR_by_state.tsv")
state_rr <- fread(fn_rr) %>%
  select(x, y, euclid_dist, min_cbsa_dist, RR, N_pairs, N_x, N_y, N_total) %>%
  rename(RR_seq = RR)

# Load preprocessed travel data — both air columns retained; models run separately per source
df_travel <- fread("data/travel_vars.tsv") %>%
  inner_join(state_rr, by = c("x", "y")) %>%
  mutate(same_state = x == y)

# Comparison of NHTS and Safegraph RRs
ggplot(df_travel) +
  aes(x = RR_trips, y = RR_move) +
  geom_point(alpha = 0.4, aes(color = same_state)) +
  geom_abline(slope = 1, color = "red", linetype = "dashed") +
  scale_x_continuous(transform = "log10", name = expression(RR["NHTS"]), limits = c(1E-3, 1E3)) +
  scale_y_continuous(transform = "log10", name = expression(RR["Safegraph"]), limits = c(1E-3, 1E3)) +
  theme_bw()
  
ggsave("figs/travel_RR_comparison.jpg", width = 6, height = 5, dpi = 300)

num_states <- nrow(df_travel) %>% sqrt()
k <- 1:(num_states^2)
i <- (k - 1) %% (num_states) + 1
j <- floor(k / num_states) + 1
k_sym <- num_states * (i - 1) + j

# Calculate correlations between sequence RR and mobility metrics
cor_tests <- list(
  NHTS     = cor.test(df_travel$RR_seq, df_travel$RR_trips, method = "spearman"),
  Air_T100 = cor.test(df_travel$RR_seq, df_travel$RR_air_t100, method = "spearman"),
  Air_DB1B = cor.test(df_travel$RR_seq, df_travel$RR_air_db1b, method = "spearman"),
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
  trips    = min(df_travel$RR_trips[df_travel$RR_trips > 0],       na.rm = TRUE),
  air_t100 = min(df_travel$RR_air_t100[df_travel$RR_air_t100 > 0], na.rm = TRUE),
  air_db1b = min(df_travel$RR_air_db1b[df_travel$RR_air_db1b > 0], na.rm = TRUE),
  move     = min(df_travel$RR_move[df_travel$RR_move > 0],          na.rm = TRUE)
)

# Normalize travel RR variables (z-score normalization)
# Distance variables (min_cbsa_dist, euclid_dist) are left unnormalized
normalize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Add normalized travel RR variables to the main df_travel dataframe
df_travel <- df_travel %>%
  mutate(
    zlog_RR_trips    = normalize(log10(RR_trips    + zero_corrections$trips)),
    zlog_RR_air_t100 = normalize(log10(RR_air_t100 + zero_corrections$air_t100)),
    zlog_RR_air_db1b = normalize(log10(RR_air_db1b + zero_corrections$air_db1b)),
    zlog_RR_move     = normalize(log10(RR_move      + zero_corrections$move))
  )

# Add small constant to zeros before taking logs and normalize travel RRs
df_model <- df_travel %>%
  mutate(
    log_RR_seq = log10(RR_seq)
  ) %>%
  filter(!is.infinite(log_RR_seq))

# Fit regression models — quick overview across both air sources
models <- list(
  NHTS          = lm(log_RR_seq ~ zlog_RR_trips,    data = df_model),
  Air_T100      = lm(log_RR_seq ~ zlog_RR_air_t100, data = df_model),
  Air_DB1B      = lm(log_RR_seq ~ zlog_RR_air_db1b, data = df_model),
  Mobility      = lm(log_RR_seq ~ zlog_RR_move,     data = df_model),
  Distance      = lm(log_RR_seq ~ min_cbsa_dist,    data = df_model),
  Combined_T100 = lm(log_RR_seq ~ zlog_RR_trips + zlog_RR_air_t100 +
                       zlog_RR_move + min_cbsa_dist, data = df_model),
  Combined_DB1B = lm(log_RR_seq ~ zlog_RR_trips + zlog_RR_air_db1b +
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

# Stratify air travel by intra- vs inter-regional (BEA regions)
REGION_DATA <- fread("data/us_states_regions.csv")

df_travel_dist <- df_travel %>%
  filter(x != y) %>%
  left_join(REGION_DATA %>% select(state, bea_reg), by = c("x" = "state")) %>%
  rename(bea_reg_x = bea_reg) %>%
  left_join(REGION_DATA %>% select(state, bea_reg), by = c("y" = "state")) %>%
  rename(bea_reg_y = bea_reg) %>%
  filter(!is.na(bea_reg_x), !is.na(bea_reg_y)) %>%
  mutate(flight_type = factor(
    ifelse(bea_reg_x == bea_reg_y, "Intra-regional", "Inter-regional"),
    levels = c("Intra-regional", "Inter-regional")
  ))

travel_comps <- list(c("Intra-regional", "Inter-regional"))

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
    tryCatch({
      boot_m <- car::Boot(m, R = 200)
      print(summary(boot_m))
      print(confint(boot_m))
    }, error = function(e) {
      cat("Bootstrap failed (likely interaction terms): ", conditionMessage(e), "\n")
    })
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

  # Calculate RMSE on log10 scale (avoids bias from high-RR outliers)
  rmse_val <- yardstick::rmse_vec(
    truth    = log10(df$RR_obs),
    estimate = log10(df$RR_hat)
  )
  rmse_val <- as.numeric(rmse_val)
  message(sprintf("log10 RMSE: %.3f", rmse_val))

  # Set plot range
  range_vals <- range(c(df$RR_obs, df$RR_hat), na.rm = TRUE)

  # Create plot
  p <- ggplot(data = df) +
    aes(x = .data$RR_hat, y = .data$RR_obs) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(transform = "log10", name = "Predicted") +
    scale_y_continuous(transform = "log10", name = "Observed") +
    ggtitle(sprintf("%s - log10 RMSE: %.3f", model_name, rmse_val)) +
    coord_equal(xlim = range_vals, ylim = range_vals) +
    theme_bw()

  return(p)
}


# Pre-assign spline plots to spacers; overwritten below if spline_flag = TRUE
slm_trips_plot       <- patchwork::plot_spacer()
slm_move_plot        <- patchwork::plot_spacer()
slm_air_t100_plot    <- patchwork::plot_spacer()
slm_air_int_t100_plot <- patchwork::plot_spacer()
slm_air_db1b_plot    <- patchwork::plot_spacer()
slm_air_int_db1b_plot <- patchwork::plot_spacer()
slm_euclid_plot      <- patchwork::plot_spacer()
slm_cbsa_plot        <- patchwork::plot_spacer()
slm_combo_plot       <- patchwork::plot_spacer()

df_model_data <- df_travel_dist %>%
  filter(x != y) %>%
  na.omit()  # Remove all rows with any NA values for car::Boot compatibility
y_obs <- df_model_data$RR_seq

# ============================================================================
# PREDICTOR CORRELATION HEATMAP
# ============================================================================
source("scripts/color_schemes.R")

pred_vars <- c(
  "zlog_RR_trips", "zlog_RR_move", "zlog_RR_air_t100", "zlog_RR_air_db1b", "min_cbsa_dist"
)
pred_labels <- c("NHTS", "SafeGraph", "T100", "DB1B", "CBSA dist.")

cor_mat <- cor(df_model_data[, ..pred_vars], use = "complete.obs", method = "pearson")
rownames(cor_mat) <- pred_labels
colnames(cor_mat) <- pred_labels

df_cor <- as.data.frame(as.table(cor_mat)) %>%
  rename(x = Var1, y = Var2, r = Freq) %>%
  mutate(
    x = factor(x, levels = pred_labels),
    y = factor(y, levels = rev(pred_labels)),
    label = sprintf("%.2f", r)
  )

p_cor <- ggplot(df_cor, aes(x = x, y = y, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.5, color = "black") +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, limits = c(-1, 1), name = "Pearson's r"
  ) +
  coord_fixed() +
  labs(
    title = "Predictor correlations",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.key.height = unit(1.2, "cm")
  )

dir.create(paste0("figs/", scenario, "/dist"), recursive = TRUE, showWarnings = FALSE)
ggsave(
  paste0("figs/", scenario, "/dist/predictor_cor_heatmap.jpg"),
  plot = p_cor, width = 5, height = 4.5, units = "in", dpi = 300
)
fn_supp_pdf <- "manuscript/figures/supp/predictor_cor_heatmap.pdf"
dir.create(dirname(fn_supp_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(fn_supp_pdf, plot = p_cor, width = 5, height = 4.5, units = "in")
message("Predictor correlation heatmap saved.")

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

if (spline_flag) {
  slm_trips <- lm(
    data = df_model_data,
    formula = log10(RR_seq) ~ bs(log10(RR_trips + zero_corrections$trips), df = 6)
  )
  slm_trips %>% model_summary(model_name = "NHTS Trips - Spline", output_file = results_file)
  slm_trips_plot <- plot_pred_RR(slm_trips, y_obs, "NHTS Trips - Spline")
}

lm_move <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_move + zero_corrections$move)
)
lm_move %>% model_summary(model_name = "Safegraph Mobility", output_file = results_file)
lm_move_plot <- plot_pred_RR(lm_move, y_obs, "Safegraph Mobility")

if (spline_flag) {
  slm_move <- lm(
    data = df_model_data,
    formula = log10(RR_seq) ~ bs(log10(RR_move + zero_corrections$move), df = 6)
  )
  slm_move %>% model_summary(calc_bootstrap = FALSE, model_name = "Safegraph Mobility - Spline", output_file = results_file)
  slm_move_plot <- plot_pred_RR(slm_move, y_obs, "Safegraph Mobility - Spline")
}

lm_air_t100 <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_air_t100 + zero_corrections$air_t100)
)
lm_air_t100 %>% model_summary(model_name = "Air Travel (T100)", output_file = results_file)
lm_air_t100_plot <- plot_pred_RR(lm_air_t100, y_obs, "Air Travel (T100)")

if (spline_flag) {
  slm_air_t100 <- lm(
    data = df_model_data,
    formula = log10(RR_seq) ~ bs(log10(RR_air_t100 + zero_corrections$air_t100), df = 6)
  )
  slm_air_t100 %>% model_summary(calc_bootstrap = FALSE, model_name = "Air Travel (T100) - Spline", output_file = results_file)
  slm_air_t100_plot <- plot_pred_RR(slm_air_t100, y_obs, "Air Travel (T100) - Spline")
}

lm_air_int_t100 <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_air_t100 + zero_corrections$air_t100):flight_type
)
lm_air_int_t100 %>% model_summary(model_name = "Air Travel (T100) : Flight Type", output_file = results_file)
lm_air_int_t100_plot <- plot_pred_RR(lm_air_int_t100, y_obs, "Air Travel (T100) : Flight Type")

if (spline_flag) {
  slm_air_int_t100 <- lm(
    data = df_model_data,
    log10(RR_seq) ~ bs(log10(RR_air_t100 + zero_corrections$air_t100), df = 6) * bs(min_cbsa_dist, df = 6)
  )
  slm_air_int_t100 %>% model_summary(calc_bootstrap = FALSE, model_name = "Air Travel (T100) x CBSA Distance - Spline", output_file = results_file)
  slm_air_int_t100_plot <- plot_pred_RR(slm_air_int_t100, y_obs, "Air Travel (T100) x CBSA Distance - Spline")
}

lm_air_db1b <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_air_db1b + zero_corrections$air_db1b)
)
lm_air_db1b %>% model_summary(model_name = "Air Travel (DB1B)", output_file = results_file)
lm_air_db1b_plot <- plot_pred_RR(lm_air_db1b, y_obs, "Air Travel (DB1B)")

if (spline_flag) {
  slm_air_db1b <- lm(
    data = df_model_data,
    formula = log10(RR_seq) ~ bs(log10(RR_air_db1b + zero_corrections$air_db1b), df = 6)
  )
  slm_air_db1b %>% model_summary(calc_bootstrap = FALSE, model_name = "Air Travel (DB1B) - Spline", output_file = results_file)
  slm_air_db1b_plot <- plot_pred_RR(slm_air_db1b, y_obs, "Air Travel (DB1B) - Spline")
}

lm_air_int_db1b <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ log10(RR_air_db1b + zero_corrections$air_db1b):flight_type
)
lm_air_int_db1b %>% model_summary(model_name = "Air Travel (DB1B) : Flight Type", output_file = results_file)
lm_air_int_db1b_plot <- plot_pred_RR(lm_air_int_db1b, y_obs, "Air Travel (DB1B) : Flight Type")

if (spline_flag) {
  slm_air_int_db1b <- lm(
    data = df_model_data,
    log10(RR_seq) ~ bs(log10(RR_air_db1b + zero_corrections$air_db1b), df = 6) * bs(min_cbsa_dist, df = 6)
  )
  slm_air_int_db1b %>% model_summary(calc_bootstrap = FALSE, model_name = "Air Travel (DB1B) x CBSA Distance - Spline", output_file = results_file)
  slm_air_int_db1b_plot <- plot_pred_RR(slm_air_int_db1b, y_obs, "Air Travel (DB1B) x CBSA Distance - Spline")
}

lm_euclid <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ euclid_dist
)
lm_euclid %>% model_summary(model_name = "Euclidean Distance", output_file = results_file)
lm_euclid_plot <- plot_pred_RR(lm_euclid, y_obs, "Euclidean Distance")

if (spline_flag) {
  slm_euclid <- lm(
    data = df_model_data,
    formula = log10(RR_seq) ~ bs(euclid_dist, df = 6)
  )
  slm_euclid %>% model_summary(calc_bootstrap = FALSE, model_name = "Euclidean Distance - Spline", output_file = results_file)
  slm_euclid_plot <- plot_pred_RR(slm_euclid, y_obs, "Euclidean Distance - Spline")
}

lm_cbsa <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ min_cbsa_dist
)
lm_cbsa %>% model_summary(model_name = "CBSA Distance", output_file = results_file)
lm_cbsa_plot <- plot_pred_RR(lm_cbsa, y_obs, "CBSA Distance")

if (spline_flag) {
  slm_cbsa <- lm(
    data = df_model_data,
    formula = log10(RR_seq) ~ bs(min_cbsa_dist, df = 6)
  )
  slm_cbsa %>% model_summary(calc_bootstrap = FALSE, model_name = "CBSA Distance - Spline", output_file = results_file)
  slm_cbsa_plot <- plot_pred_RR(slm_cbsa, y_obs, "CBSA Distance - Spline")
}

# Combined Regression — DB1B removed (suppressor); T100 retained
lm_combo <- lm(
  data = df_model_data,
  formula = log10(RR_seq) ~ zlog_RR_trips +
    zlog_RR_move +
    min_cbsa_dist +
    zlog_RR_air_t100 +
    zlog_RR_air_t100:flight_type
)
lm_combo %>% model_summary(model_name = "Multivariate Model (normalized)", output_file = results_file)
lm_combo_plot <- plot_pred_RR(lm_combo, y_obs, "Multivariate Model (normalized)")

# ============================================================================
# STEPWISE BIC VARIABLE SELECTION
# ============================================================================
message("\n", strrep("=", 80))
message("STEPWISE BIC VARIABLE SELECTION")
message(strrep("=", 80), "\n")

lm_bic <- MASS::stepAIC(lm_combo, direction = "both", k = log(nrow(df_model_data)), trace = TRUE)

sink(results_file, append = TRUE)
cat(paste0("\n", strrep("=", 80), "\n"))
cat("STEPWISE BIC SELECTED MODEL\n")
cat(paste0(strrep("=", 80), "\n\n"))
cat("Selected formula:\n")
print(formula(lm_bic))
cat("\n")
print(summary(lm_bic))
sink()

message("\nBIC-selected formula:")
print(formula(lm_bic))

lm_bic %>% model_summary(model_name = "BIC-Selected Model", output_file = results_file)
lm_bic_plot <- plot_pred_RR(lm_bic, y_obs, "BIC-Selected Model")

# ============================================================================
# VARIANCE PARTITIONING — run_varpart helper + calls for all/intra/inter
# ============================================================================
suppressPackageStartupMessages({
  library(relaimpo)
  library(vegan)
})

GROUP_PATTERNS <- list(
  "Ground travel\n(NHTS)"  = "RR_trips",
  "Mobility\n(SafeGraph)"  = "RR_move",
  "Distance\n(CBSA)"       = "cbsa_dist",
  "Air travel\n(T100)"     = "air_t100"
)
GROUP_COLORS <- c("steelblue", "seagreen", "firebrick", "darkorange")

run_varpart <- function(lm_fit, label, out_prefix) {
  message(sprintf("\n--- Variance partitioning: %s ---", label))
  mm <- model.matrix(lm_fit)[, -1, drop = FALSE]
  cn <- colnames(mm)
  Y  <- model.response(model.frame(lm_fit))

  # Detect which groups are represented in this model's matrix
  X_list   <- Filter(function(x) ncol(x) > 0,
                     lapply(GROUP_PATTERNS, function(p) mm[, grepl(p, cn), drop = FALSE]))
  grp_cols <- unlist(lapply(GROUP_PATTERNS, function(p) sum(grepl(p, cn))))
  present  <- names(grp_cols)[grp_cols > 0]
  colors   <- GROUP_COLORS[match(present, names(GROUP_PATTERNS))]

  # LMG via expanded plain model matrix
  df_mm   <- as.data.frame(mm); df_mm$y <- Y
  lm_exp  <- lm(y ~ ., data = df_mm)
  lmg_res <- calc.relimp(lm_exp, type = "lmg", rela = FALSE)

  df_lmg <- data.frame(term = names(lmg_res@lmg), lmg = as.numeric(lmg_res@lmg)) %>%
    mutate(group = case_when(
      grepl("RR_trips",  term) ~ "Ground travel\n(NHTS)",
      grepl("RR_move",   term) ~ "Mobility\n(SafeGraph)",
      grepl("cbsa_dist", term) ~ "Distance\n(CBSA)",
      grepl("air_t100",  term) ~ "Air travel\n(T100)"
    )) %>%
    group_by(group) %>%
    summarise(lmg = sum(lmg), .groups = "drop") %>%
    arrange(desc(lmg))

  message("LMG by group:"); print(df_lmg)
  message(sprintf("Total R² = %.4f", sum(df_lmg$lmg)))

  p_lmg <- ggplot(df_lmg, aes(x = reorder(group, lmg), y = lmg)) +
    geom_col(fill = "steelblue", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", lmg * 100)), hjust = -0.15, size = 4.5) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, max(df_lmg$lmg) * 1.25),
                       name = expression("Relative importance (LMG " * R^2 * ")")) +
    coord_flip() +
    labs(title = sprintf("Relative importance — %s", label),
         subtitle = sprintf("BIC-selected model  |  Total R² = %.3f", sum(df_lmg$lmg)),
         x = NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(color = "grey40"))

  ggsave(paste0(out_prefix, "_lmg.jpg"), plot = p_lmg,
         width = 7, height = 4, dpi = 300)

  # Variance partitioning Venn
  if (length(X_list) >= 2) {
    short_names <- c(
      "Ground travel\n(NHTS)" = "NHTS",
      "Mobility\n(SafeGraph)" = "SafeGraph",
      "Distance\n(CBSA)"      = "CBSA dist.",
      "Air travel\n(T100)"    = "T100"
    )
    vp_names <- unname(short_names[names(X_list)])

    vp <- do.call(varpart, c(list(Y), unname(X_list)))
    message("\nVariance partitioning fractions:"); print(vp)
    jpeg(paste0(out_prefix, "_varpart.jpg"), width = 10, height = 9, units = "in", res = 300)
    plot(vp, digits = 2, Xnames = vp_names, bg = colors, alpha = 180, cex = 1.1)
    title(sprintf("Variance partitioning: %s", label), cex.main = 1.3)
    dev.off()
    message(sprintf("Plots saved to %s_lmg.jpg / %s_varpart.jpg", out_prefix, out_prefix))
  }
}

message("\n", strrep("=", 80))
message("VARIANCE PARTITIONING (BIC-selected model)")
message(strrep("=", 80))
run_varpart(lm_bic,
            label      = "All pairs",
            out_prefix = paste0("figs/", scenario, "/dist/all"))

# ============================================================================
# STRATIFIED MODELS: INTRA- vs INTER-REGIONAL
# ============================================================================
message("\n", strrep("=", 80))
message("STRATIFIED MODELS BY FLIGHT TYPE")
message(strrep("=", 80), "\n")

# Within each stratum flight_type is constant, so interactions collapse to
# simple main effects — formula is the same for both strata.
# Use subset= so stepAIC can find df_model_data in global scope when refitting.
formula_strat <- log10(RR_seq) ~ zlog_RR_trips + zlog_RR_move + min_cbsa_dist +
  zlog_RR_air_t100

strat_results <- list()

for (stratum in c("Intra-regional", "Inter-regional")) {
  n_s <- sum(df_model_data$flight_type == stratum)
  message(sprintf("\n--- %s  (n = %d) ---", stratum, n_s))

  lm_s <- lm(formula_strat, data = df_model_data,
              subset = flight_type == stratum)

  message("Full model:")
  print(summary(lm_s)$coefficients)
  message(sprintf("Adj. R² = %.3f", summary(lm_s)$adj.r.squared))

  message("\nBIC stepwise:")
  lm_s_bic <- MASS::stepAIC(lm_s, direction = "both",
                              k = log(n_s), trace = TRUE)

  message("\nBIC-selected formula:")
  print(formula(lm_s_bic))
  message(sprintf("Adj. R² (BIC model) = %.3f", summary(lm_s_bic)$adj.r.squared))

  robust_s <- lmtest::coeftest(lm_s_bic, vcov = sandwich::vcovHC(lm_s_bic, type = "HC1"))
  robust_ci_s <- lmtest::coefci(lm_s_bic, vcov = sandwich::vcovHC(lm_s_bic, type = "HC1"))

  sink(results_file, append = TRUE)
  cat(paste0("\n", strrep("=", 80), "\n"))
  cat(paste0("STRATIFIED MODEL: ", stratum, "  (n = ", n_s, ")\n"))
  cat(paste0(strrep("=", 80), "\n\n"))
  cat("--- Full model ---\n")
  print(summary(lm_s))
  cat("\n--- BIC-selected model ---\n")
  cat("Formula: "); print(formula(lm_s_bic))
  print(summary(lm_s_bic))
  if (length(attr(terms(lm_s_bic), "term.labels")) >= 2) {
    cat("\nVARIANCE INFLATION FACTORS (VIF)\n")
    print(car::vif(lm_s_bic))
  }
  cat("\nROBUST COEFFICIENT TEST (HC1 ROBUST SEs)\n")
  print(robust_s)
  cat("\n95% Confidence Intervals (Robust SEs):\n")
  print(robust_ci_s)
  sink()

  run_varpart(lm_s_bic,
              label      = stratum,
              out_prefix = paste0("figs/", scenario, "/dist/",
                                  tolower(gsub("-| ", "_", stratum))))

  strat_results[[stratum]] <- list(stratum = stratum, n = n_s,
                                    lm_full = lm_s, lm_bic = lm_s_bic)
}

# Coefficient comparison table across strata
message("\n--- Coefficient comparison across strata (BIC-selected models) ---")
coef_table <- map_dfr(strat_results, function(res) {
  cf <- summary(res$lm_bic)$coefficients
  as.data.frame(cf) %>%
    rownames_to_column("term") %>%
    mutate(stratum = res$stratum, adj_r2 = summary(res$lm_bic)$adj.r.squared)
})
print(coef_table)

# ============================================================================
# NEGATIVE BINOMIAL GLMs + TENSOR PRODUCT GAMs
# ============================================================================
source("scripts/state_tensor_plots.R")

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
  dplyr::select(x, y, RR_seq, predicted_RR, residual, studentized_residual,
         cooks_distance, RR_trips, RR_air_t100, RR_air_db1b, RR_move, min_cbsa_dist) %>%
  arrange(desc(studentized_residual))

message("\nTop enriched state pairs (observed > expected):")
print(head(df_enriched, 20))

# Top depleted pairs (transmission lower than expected from travel)
df_depleted <- df_anomalies %>%
  filter(anomaly_type == "Depleted") %>%
  dplyr::select(x, y, RR_seq, predicted_RR, residual, studentized_residual,
         cooks_distance, RR_trips, RR_air_t100, RR_air_db1b, RR_move, min_cbsa_dist) %>%
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
       plot = plot_residuals, width = 10, height = 7, dpi = 300)

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
       plot = plot_cooks, width = 10, height = 6, dpi = 300)

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
       plot = plot_qq, width = 8, height = 6, dpi = 300)

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
       plot = plot_obs_pred_anomaly, width = 10, height = 8, dpi = 300)

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
    min_cbsa_dist +
    zlog_RR_air_t100 +
    zlog_RR_air_t100:flight_type
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
       plot = comparison_plot, width = 14, height = 6, dpi = 300)

message("\nComparison plot saved to figs/", scenario, "/model_comparison_with_without_anomalies.jpg")

# Spline combined regression
if (spline_flag) {
  slm_combo <- lm(
    data = df_model_data,
    formula = log10(RR_seq) ~ zlog_RR_trips +
      zlog_RR_move +
      min_cbsa_dist +
      bs(zlog_RR_air_t100, df = 6) * bs(min_cbsa_dist, df = 6)
  )
  slm_combo %>% model_summary(calc_bootstrap = FALSE, model_name = "Multivariate Model - Spline (normalized)", output_file = results_file)
  slm_combo_plot <- plot_pred_RR(slm_combo, y_obs, "Multivariate Model - Spline (normalized)")
}

if (spline_flag) {
  OE_plots <- patchwork::wrap_plots(
    lm_trips_plot,    slm_trips_plot,    lm_move_plot,         slm_move_plot,
    lm_air_t100_plot, slm_air_t100_plot, lm_air_int_t100_plot, slm_air_int_t100_plot,
    lm_air_db1b_plot, slm_air_db1b_plot, lm_air_int_db1b_plot, slm_air_int_db1b_plot,
    lm_euclid_plot,   slm_euclid_plot,   lm_cbsa_plot,         slm_cbsa_plot,
    lm_combo_plot,    slm_combo_plot,    lm_bic_plot,          patchwork::plot_spacer(),
    ncol = 4, nrow = 5
  )
  oe_w <- 22; oe_h <- 25
} else {
  OE_plots <- patchwork::wrap_plots(
    lm_trips_plot,    lm_move_plot,
    lm_air_t100_plot, lm_air_int_t100_plot,
    lm_air_db1b_plot, lm_air_int_db1b_plot,
    lm_euclid_plot,   lm_cbsa_plot,
    lm_combo_plot,    lm_bic_plot,
    ncol = 2, nrow = 5
  )
  oe_w <- 11; oe_h <- 25
}
ggsave(paste0("figs/", scenario, "/dist/regression_fit.jpg"),
  plot = OE_plots, width = oe_w, height = oe_h, units = "in", dpi = 300
)
