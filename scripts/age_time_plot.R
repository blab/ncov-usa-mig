library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = 'CAM_1000', help = 'Scenario name')
  return(parser$parse_args())
}

# Function to generate all three plots for a given dataset
generate_age_time_plots <- function(df_data, fig_path, title_suffix = "") {
  # Define age group order
  age_order <- c("0-5y", "6-11y", "12-17y", "18-25y", "26-45y", "46-65y", "65-80y", "81y+")

  # Set factor levels
  df_data <- df_data %>%
    mutate(
      x = factor(x, levels = age_order),
      y = factor(y, levels = age_order)
    )

  # Plot 1: RR over time - include diagonal (all pairs)
  p1 <- ggplot(df_data, aes(x = date, y = RR, color = y, group = y)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black", alpha = 0.5) +
    facet_wrap(~ x, scales = "free_y", nrow = 2, ncol = 4) +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
    scale_y_log10() +
    labs(
      title = paste0("Relative Risk (RR) Over Time by Age Group", title_suffix),
      x = "Date",
      y = "RR (log scale)",
      color = "Paired Age Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "bottom"
    )

  ggsave(
    filename = paste0(fig_path, "age_RR_time_faceted.png"),
    plot = p1,
    width = 16,
    height = 8,
    dpi = 192
  )

  # Plot 2: nRR over time - exclude diagonal
  plot_data_nrr <- df_data %>%
    filter(x != y)

  p2 <- ggplot(plot_data_nrr, aes(x = date, y = nRR, color = y, group = y)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
    facet_wrap(~ x, scales = "free_y", nrow = 2, ncol = 4) +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
    labs(
      title = paste0("Diagonal Normalized Relative Risk (nRR) Over Time by Age Group", title_suffix),
      x = "Date",
      y = "nRR",
      color = "Paired Age Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "bottom"
    )

  ggsave(
    filename = paste0(fig_path, "age_nRR_time_faceted.png"),
    plot = p2,
    width = 16,
    height = 8,
    dpi = 192
  )

  # Plot 3: nRR_fixed over time - include diagonal (all pairs)
  p3 <- ggplot(df_data, aes(x = date, y = nRR_fixed, color = y, group = y)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "black", alpha = 0.5) +
    facet_wrap(~ x, scales = "free_y", nrow = 2, ncol = 4) +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
    labs(
      title = paste0("Baseline Normalized RR Over Time by Age Group (vs 26-45y)", title_suffix),
      x = "Date",
      y = "nRR (normalized to 26-45y)",
      color = "Paired Age Group"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "bottom"
    )

  ggsave(
    filename = paste0(fig_path, "age_nRR_fixed_time_faceted_all.png"),
    plot = p3,
    width = 16,
    height = 8,
    dpi = 192
  )
}

args <- collect_args()
scenario <- args$scenario

# 1. Generate plots for all countries combined
cat("Generating plots for all countries combined...\n")
data_path_all <- paste0("results/", scenario, "/time_age/df_RR_by_time_age_series.tsv")
fig_path_all <- paste0("figs/", scenario, "/age_time/all_countries/")
dir.create(fig_path_all, recursive = TRUE, showWarnings = FALSE)

df_age_time_all <- read_tsv(data_path_all)
generate_age_time_plots(df_age_time_all, fig_path_all, title_suffix = " (All Countries)")

cat("Plots saved to:", fig_path_all, "\n")

# 2. Generate plots for each country
cat("Generating plots for individual countries...\n")
data_path_countries <- paste0("results/", scenario, "/time_age/df_RR_by_time_age_countries.tsv")

# Check if country-specific data exists
if (file.exists(data_path_countries)) {
  df_age_time_countries <- read_tsv(data_path_countries)

  # Get unique countries
  countries <- df_age_time_countries %>%
    distinct(country) %>%
    pull(country)

  for (ctry in countries) {
    cat("  Processing country:", ctry, "\n")

    # Filter data for this country
    df_country <- df_age_time_countries %>%
      filter(country == ctry)

    # Create country-specific figure directory
    fig_path_country <- paste0("figs/", scenario, "/age_time/", ctry, "/")
    dir.create(fig_path_country, recursive = TRUE, showWarnings = FALSE)

    # Generate plots
    generate_age_time_plots(df_country, fig_path_country, title_suffix = paste0(" (", ctry, ")"))

    cat("  Plots saved to:", fig_path_country, "\n")
  }
} else {
  cat("Country-specific data not found at:", data_path_countries, "\n")
  cat("Skipping country-specific plots.\n")
}

cat("\nAll plots completed.\n")
