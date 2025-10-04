# File: state_time_dist.R
# Author(s): Amin Bemanian
# Date: 2025-09-25
# Description: Analyze and visualize temporal patterns in state-level relative risk (RR) 
#              and its relationship with geographic distances and travel patterns
# 
# Arguments:
#   None - Uses fixed SCENARIO variable for file paths
#
# Input files:
#   - data/nb_dist_states.tsv: State neighbor distances matrix
#   - data/travel_vars.tsv: Preprocessed travel variables
#   - results/{SCENARIO}/time_state/df_state_rr_snap.tsv: Snapshot RR data
#   - results/{SCENARIO}/time_state/df_state_rr_series.tsv: Time series RR data
#
# Output:
#   Generated plots in figs/{SCENARIO}/time/:
#   - rr_boxplot_time.png: RR distribution over time
#   - nb_boxplot_time.png: RR by neighbor rank
#   - euclid_dist_time.png: RR vs geographic distance
#   - nhts_time.png: RR vs travel patterns
#
# Dependencies: tidyverse, data.table, ggplot2, tictoc, scales, maps, sf

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(tictoc)
  library(scales)
  library(maps)
  library(sf)
})

# Import specific functions to avoid namespace conflicts
`%>%` <- dplyr::`%>%`
select <- dplyr::select
left_join <- dplyr::left_join
rename <- dplyr::rename
filter <- dplyr::filter
group_by <- dplyr::group_by
summarise <- dplyr::summarise
ungroup <- dplyr::ungroup
join_by <- dplyr::join_by
across <- dplyr::across
all_of <- dplyr::all_of
pivot_longer <- tidyr::pivot_longer
pivot_wider <- tidyr::pivot_wider
starts_with <- dplyr::starts_with
mutate <- dplyr::mutate
str_remove <- stringr::str_remove
everything <- dplyr::everything
map2 <- purrr::map2
map_dbl <- purrr::map_dbl
group_modify <- dplyr::group_modify
bind_rows <- dplyr::bind_rows
select <- dplyr::select

# Prevent no visible binding notes
utils::globalVariables(c(
  "x", "y", "state_x", "state_y", "RR_trips", "nRR",
  "date", "date_median_nRR", "nb_dist", "euclid_dist",
  "RR_nhts", "RR_move", "RR_air", "variable", "correlation",
  "metric", "value", "n", "ci_lower", "ci_upper", "significant",
  "cor", "ci", ".x", ".y", "process_variable",
  # Additional variables needed for correlation functions
  "get", "correlation", "variable", "ci",
  # Variables for ggplot
  "x", "y", "nRR", "date", "correlation", "variable",
  "ci_lower", "ci_upper", "significant"
))

df_state_distances <- fread("data/nb_dist_states.tsv")
df_cbsa_distances <- fread("data/state_cbsa_dist.tsv")
df_travel <- fread("data/travel_vars.tsv")

SCENARIO <- "CAM_1000"
fig_path <- paste0("figs/",SCENARIO,"/time/")

# Add the normalization function from state_time_rr_analysis.R
normalized_state_rr <- function(df_rr){
  # First: filter to the diagonal (RR(x,x)) entries
  rr_diag <- df_rr %>%
    filter(x == y) %>%
    select(date, state = x, RR_diag = RR, N_diag = N_pairs) %>%
    mutate(RR_diag = ifelse(N_diag == 0,NA,RR_diag))
  # Join to get RR(x,x) and RR(y,y) for each row
  df_rr <- df_rr %>%
    left_join(rr_diag, by = c("date", "x" = "state")) %>%
    rename(RR_xx = RR_diag) %>%
    left_join(rr_diag, by = c("date", "y" = "state")) %>%
    rename(RR_yy = RR_diag) %>%
    mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
return(df_rr)
}

# Modified normalization function for travel data (without N_pairs)
normalized_travel_rr <- function(df_rr){
  # First: filter to the diagonal (RR(x,x)) entries
  rr_diag <- df_rr %>%
    filter(x == y) %>%
    select(date, state = x, RR_diag = RR) %>%
    mutate(RR_diag = ifelse(is.na(RR_diag), NA, RR_diag))
  
  # Join to get RR(x,x) and RR(y,y) for each row
  df_rr <- df_rr %>%
    left_join(rr_diag, by = c("date", "x" = "state")) %>%
    rename(RR_xx = RR_diag) %>%
    left_join(rr_diag, by = c("date", "y" = "state")) %>%
    rename(RR_yy = RR_diag) %>%
    mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
  return(df_rr)
}

# Function to normalize travel RR variables
normalize_travel_vars <- function(df) {
  # Normalize RR_trips
  if("RR_trips" %in% names(df)) {
    df_trips <- df %>%
      select(date, x, y, RR = RR_trips) %>%
      filter(!is.na(RR)) %>%
      normalized_travel_rr() %>%
      select(date, x, y, nRR_trips = nRR)
    df <- left_join(df, df_trips, by = c("date", "x", "y"))
  }
  
  # Normalize RR_move
  if("RR_move" %in% names(df)) {
    df_move <- df %>%
      select(date, x, y, RR = RR_move) %>%
      filter(!is.na(RR)) %>%
      normalized_travel_rr() %>%
      select(date, x, y, nRR_move = nRR)
    df <- left_join(df, df_move, by = c("date", "x", "y"))
  }
  
  return(df)
}

attach_distances <- function(df) {
  # Join neighbor distances
  df <- left_join(
    df,
    df_state_distances,
    by = c("x" = "state_x", "y" = "state_y")
  )
  
  df <- left_join(
    df,
    df_cbsa_distances,
    by = c("x", "y")
  )
  
  # Join travel data and rename
  df <- left_join(
    df,
    df_travel,
    by = c("x", "y")
  )
  
  df <- df %>%
    mutate(
       RR_air_short = ifelse(min_cbsa_dist < 300, RR_air, NA),
       RR_air_med = ifelse(min_cbsa_dist >= 300 & min_cbsa_dist < 1800, RR_air, NA),
       RR_air_long = ifelse(min_cbsa_dist >= 1800, RR_air, NA)
    ) %>%
    normalize_travel_vars()
    
  return(df)
}

state_rr_snap <- read_tsv(paste0("results/",
                                 SCENARIO,
                                 "/time_state/df_state_rr_snap.tsv")) %>%
  attach_distances()
state_rr_series <- read_tsv(paste0("results/",
                                   SCENARIO,
                                   "/time_state/df_state_rr_series.tsv")) %>%
  attach_distances()

date_medians <- state_rr_snap%>%
  filter(x != y) %>%
  group_by(date) %>%
  summarise(date_median_nRR = median(nRR, na.rm = TRUE)) %>%
  ungroup()

# Create standardized time series boxplot function
plot_rr_boxplot <- function(data, date_medians, filename) {
  p <- ggplot(data %>%
                filter(x != y) %>%
                left_join(date_medians, by = "date"),
              aes(x = date, y = nRR, group = date)) +
    geom_jitter(alpha = 0.1) +
    geom_boxplot(aes(fill = date_median_nRR), outliers = FALSE, alpha = 0.9) +
    scale_fill_gradient(low = "gold", high = "red2",
                       transform = "log10",
                       limits = c(0.0001, 0.07),
                       oob = scales::squish,
                       na.value = "gold",
                       labels = label_number(accuracy = 0.01)) +
    ylim(0, 0.5) + 
    labs(title = "Interstate Normalized RR over Time",
         y = "nRR",
         x = "Date of RR Snapshot",
         fill = "Median nRR") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(filename = paste0(fig_path, filename),
         plot = p,
         width = 7,
         height = 5,
         units = "in",
         dpi = 192,
         create.dir = TRUE)
  
  return(p)
}

# Create standardized neighbor rank boxplot function
plot_neighbor_rank <- function(data, filename) {
  p <- ggplot(data %>% filter(x != y),
              aes(x = nb_dist, y = nRR, group = interaction(nb_dist, date))) +
    geom_jitter(alpha = 0.1) +
    geom_boxplot(aes(fill = as.factor(date)),
                 alpha = 0.9,
                 outliers = FALSE,
                 show.legend = FALSE) +
    theme_bw() +
    ylim(0, 0.5) +  
    xlim(0.5, 7.5) +
    facet_wrap(vars(date), nrow = 4, dir = "v") +
    labs(title = "Neighbor Rank and Normalized RR",
         x = "Neighbor Rank (Queen Adjacency)",
         y = "Normalized RR")
  
  ggsave(filename = paste0(fig_path, filename),
         plot = p,
         width = 2400,
         height = 2400,
         units = "px",
         dpi = 192)
  
  return(p)
}

# Generate plots
plot_rr_boxplot(state_rr_snap, date_medians, "rr_boxplot_time.png")
plot_neighbor_rank(state_rr_snap, "nb_boxplot_time.png")

# Create standardized time series plot function
plot_rr_relationship <- function(data, x_var, title, x_label, filename, 
                               log_x = FALSE, x_limits = NULL) {
  
  # Filter data and check if variable exists
  plot_data <- data %>% filter(x != y)
  
  if (!x_var %in% names(plot_data)) {
    warning(paste("Variable", x_var, "not found in data. Skipping plot."))
    return(NULL)
  }
  
  # Check if there's enough data for each facet
  data_summary <- plot_data %>%
    group_by(date) %>%
    summarise(
      n_valid = sum(!is.na(.data[[x_var]]) & !is.na(nRR)),
      .groups = 'drop'
    )
  
  # Only proceed if we have enough data points
  min_points <- max(data_summary$n_valid)
  if (min_points < 10) {
    warning(paste("Not enough data points for", x_var, ". Maximum points in any date:", min_points))
    return(NULL)
  }
  
  p <- ggplot(plot_data, aes(x = .data[[x_var]], y = nRR)) +
    geom_point(alpha = 0.2,
               aes(color = as.factor(date)),
               show.legend = FALSE) +
    theme_bw() +
    facet_wrap(vars(date), nrow = 4, dir = "v") +
    ylim(0, 0.5) +
    labs(title = title,
         x = x_label,
         y = "Normalized RR")
  
  # Add smooth line only if we have enough points, with reduced complexity
  tryCatch({
    # Use fewer knots and simpler basis
    k_value <- min(5, floor(min_points / 3))  # Adaptive k based on data
    if (k_value >= 3) {
      p <- p + geom_smooth(method = "gam", 
                          formula = y ~ s(x, k = k_value, bs = "cs"),
                          se = FALSE)
    } else {
      # Fall back to loess for very sparse data
      p <- p + geom_smooth(method = "loess", se = FALSE)
    }
  }, error = function(e) {
    warning(paste("Could not add smooth line for", x_var, ":", e$message))
  })
  
  # Add log scale with limits if requested
  if (log_x) {
    p <- p + scale_x_log10(limits = x_limits)
  } else if (!is.null(x_limits)) {
    p <- p + xlim(x_limits[1], x_limits[2])
  }
  
  # Save plot
  ggsave(filename = paste0(fig_path, filename),
         plot = p,
         width = 2400,
         height = 2400,
         units = "px",
         dpi = 192)
  
  return(p)
}

# Plot Euclidean distance relationship
plot_rr_relationship(
  data = state_rr_snap,
  x_var = "min_cbsa_dist",
  title = "State Distance and Normalized RR",
  x_label = "Distance (km)",
  filename = "min_cbsa_dist_time.png",
  x_limits = c(0, 5000)
)

# Plot travel relationship plots
plot_rr_relationship(
  data = state_rr_snap,
  x_var = "RR_trips",
  title = "Travel RR and Normalized RR",
  x_label = "Travel RR",
  filename = "nhts_time.png",
  log_x = TRUE,
  x_limits = c(0.0001, 10)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "RR_move",
  title = "Mobility RR and Normalized RR",
  x_label = "Mobility RR",
  filename = "safegraph_time.png",
  log_x = TRUE,
  x_limits = c(0.01, 10)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "RR_air",
  title = "Air Travel RR and Normalized RR",
  x_label = "Air Travel RR",
  filename = "air_travel_time.png",
  log_x = TRUE,
  x_limits = c(0.01, 10)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "RR_air_short",
  title = "Air Travel RR (Short) and Normalized RR",
  x_label = "Air Travel RR (Short)",
  filename = "air_travel_short_time.png",
  log_x = TRUE,
  x_limits = c(0.01, 10)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "RR_air_med",
  title = "Air Travel RR (Medium) and Normalized RR",
  x_label = "Air Travel RR (Medium)",
  filename = "air_travel_med_time.png",
  log_x = TRUE,
  x_limits = c(0.01, 10)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "RR_air_long",
  title = "Air Travel RR (Long) and Normalized RR",
  x_label = "Air Travel RR (Long)",
  filename = "air_travel_long_time.png",
  log_x = TRUE,
  x_limits = c(0.01, 10)
)

# Add new normalized travel plots
plot_rr_relationship(
  data = state_rr_snap,
  x_var = "nRR_trips",
  title = "Normalized Travel RR and Normalized RR",
  x_label = "Normalized Travel RR",
  filename = "nhts_time_normalized.png",
  log_x = FALSE,
  x_limits = c(0, 5)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "nRR_move",
  title = "Normalized Mobility RR and Normalized RR",
  x_label = "Normalized Mobility RR",
  filename = "safegraph_time_normalized.png",
  log_x = FALSE,
  x_limits = c(0, 5)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "nRR_air",
  title = "Normalized Air Travel RR and Normalized RR",
  x_label = "Normalized Air Travel RR",
  filename = "air_travel_time_normalized.png",
  log_x = FALSE,
  x_limits = c(0, 5)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "nRR_air_short",
  title = "Normalized Air Travel RR (Short) and Normalized RR",
  x_label = "Normalized Air Travel RR (Short)",
  filename = "air_travel_short_time_normalized.png",
  log_x = FALSE,
  x_limits = c(0, 5)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "nRR_air_med",
  title = "Normalized Air Travel RR (Medium) and Normalized RR",
  x_label = "Normalized Air Travel RR (Medium)",
  filename = "air_travel_med_time_normalized.png",
  log_x = FALSE,
  x_limits = c(0, 5)
)

plot_rr_relationship(
  data = state_rr_snap,
  x_var = "nRR_air_long",
  title = "Normalized Air Travel RR (Long) and Normalized RR",
  x_label = "Normalized Air Travel RR (Long)",
  filename = "air_travel_long_time_normalized.png",
  log_x = FALSE,
  x_limits = c(0, 5)
)

# Function to calculate correlation CI using Fisher z-transformation
calculate_correlation_ci <- function(r, n, conf.level = 0.95) {
  # Fisher's z-transformation
  z <- 0.5 * log((1 + r) / (1 - r))
  
  # Standard error of z
  se <- 1 / sqrt(n - 3)
  
  # Critical value
  crit <- qnorm((1 + conf.level) / 2)
  
  # Confidence interval for z
  ci_lower_z <- z - crit * se
  ci_upper_z <- z + crit * se
  
  # Transform back to correlation scale
  ci_lower_r <- (exp(2 * ci_lower_z) - 1) / (exp(2 * ci_lower_z) + 1)
  ci_upper_r <- (exp(2 * ci_upper_z) - 1) / (exp(2 * ci_upper_z) + 1)
  
  return(c(ci_lower = ci_lower_r, ci_upper = ci_upper_r))
}

# Create function to calculate and plot time series correlations
plot_correlation_series <- function(data, var_list, var_labels = NULL, filename,
                                  conf.level = 0.95) {
  # Initialize empty list for results
  results_list <- vector("list", length(var_list))
  names(results_list) <- var_list
  
  # Process each variable
  for (i in seq_along(var_list)) {
    var <- var_list[i]
    var_data <- data %>% filter(x != y)
    
    # Calculate correlations for each date
    dates <- unique(var_data$date)
    corr_results <- data.frame(
      date = dates,
      correlation = numeric(length(dates)),
      n = numeric(length(dates)),
      ci_lower = numeric(length(dates)),
      ci_upper = numeric(length(dates)),
      significant = logical(length(dates)),
      variable = var
    )
    
    # Calculate correlations and CIs for each date
    for (j in seq_along(dates)) {
      date_data <- var_data[var_data$date == dates[j], ]
      valid_idx <- !is.na(date_data[[var]]) & !is.na(date_data$nRR)
      
      if (sum(valid_idx) >= 3) {
        corr_results$correlation[j] <- cor(
          date_data[[var]][valid_idx],
          date_data$nRR[valid_idx],
          method = "spearman"
        )
        corr_results$n[j] <- sum(valid_idx)
        ci <- calculate_correlation_ci(
          corr_results$correlation[j],
          corr_results$n[j],
          conf.level
        )
        corr_results$ci_lower[j] <- ci[1]
        corr_results$ci_upper[j] <- ci[2]
        corr_results$significant[j] <- (ci[1] > 0) | (ci[2] < 0)
      } else {
        corr_results$correlation[j] <- NA
        corr_results$n[j] <- sum(valid_idx)
        corr_results$ci_lower[j] <- NA
        corr_results$ci_upper[j] <- NA
        corr_results$significant[j] <- FALSE
      }
    }
    
    results_list[[i]] <- corr_results
  }
  
  # Combine results
  results_df <- do.call(rbind, results_list)
  
  # Add variable labels if provided
  if (!is.null(var_labels)) {
    results_df$variable <- factor(results_df$variable,
                                levels = var_list,
                                labels = var_labels)
  } else {
    results_df$variable <- factor(results_df$variable,
                                levels = var_list)
  }
  
  # Create plot
  p <- ggplot(results_df, aes(x = date, y = correlation, color = variable)) +
    # Add confidence interval ribbons
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper, fill = variable),
      alpha = 0.2,
      color = NA
    ) +
    # Add lines and points
    geom_line(linewidth = 1) +
    geom_point(aes(size = significant), alpha = 0.8) +
    # Add horizontal line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    # Customize scales
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.5)) +
    guides(size = "none") +  # Hide size legend
    # Theme and labels
    theme_bw() +
    labs(
      title = "Spearman Correlation with Normalized RR Over Time",
      subtitle = paste0(
        conf.level * 100,
        "% Confidence Intervals (Larger points indicate significant correlations)"
      ),
      x = "Date",
      y = "Spearman Correlation",
      color = "Variable",
      fill = "Variable"
    ) +
    ylim(-1, 1) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
  
  # Save plot
  ggsave(filename = paste0(fig_path, filename),
         plot = p,
         width = 10,
         height = 6,
         units = "in",
         dpi = 192)
  
  return(p)
  
  # Create plot with confidence intervals and significance indicators
  p <- ggplot(results_df, aes(x = date, y = correlation, color = variable)) +
    # Add confidence interval ribbons
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper, fill = variable),
      alpha = 0.2,
      color = NA
    ) +
    # Add lines and points
    geom_line(linewidth = 1) +
    geom_point(aes(size = significant), alpha = 0.8) +
    # Add horizontal line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    # Customize scales
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.5)) +
    guides(size = "none") +  # Hide size legend
    # Theme and labels
    theme_bw() +
    labs(
      title = "Spearman Correlation with Normalized RR Over Time",
      subtitle = paste0(
        conf.level * 100,
        "% Confidence Intervals (Larger points indicate significant correlations)"
      ),
      x = "Date",
      y = "Spearman Correlation",
      color = "Variable",
      fill = "Variable"
    ) +
    ylim(-1, 1) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
  
  # Save plot
  ggsave(filename = paste0(fig_path, filename),
         plot = p,
         width = 10,
         height = 6,
         units = "in",
         dpi = 192)
  
  return(p)
}

# Plot correlation time series
plot_correlation_series(
  data = state_rr_series,
  var_list = c("min_cbsa_dist", "RR_trips", "RR_move", "RR_air"),
  var_labels = c("Geographic Distance", "Travel", "Mobility", "Air Travel"),
  filename = "correlations_time_series.png"
)

plot_correlation_series(
  data = state_rr_series,
  var_list = c("min_cbsa_dist", "RR_trips", "RR_move", "RR_air_short", "RR_air_med", "RR_air_long"),
  var_labels = c("Geographic Distance", "Travel", "Mobility", "Air Travel (Short)", "Air Travel (Medium)", "Air Travel (Long)"),
  filename = "correlations_time_series_split_air.png"
)

plot_correlation_series(
  data = state_rr_series,
  var_list = c("min_cbsa_dist", "nRR_trips", "nRR_move"),
  var_labels = c("Geographic Distance", "Normalized Travel", "Normalized Mobility"),
  filename = "correlations_time_series_normalized.png"
)
