library(tidyverse)
library(argparse)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = 'CAM_1000', help = 'Scenario name')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# Define paths
data_path <- paste0("results/", scenario, "/df_RR_by_age_time_series.tsv")
fig_path <- paste0("figs/", scenario, "/age_time/")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

# Read data
df_age_time <- read_tsv(data_path)

# Calculate normalized RR (nRR)
# nRR = RR(x,y) / mean(RR(x,x), RR(y,y))
normalized_age_rr <- function(df_rr){
  # First: filter to the diagonal (RR(x,x)) entries
  rr_diag <- df_rr %>%
    filter(x == y) %>%
    select(date, age = x, RR_diag = RR, N_diag = N_pairs) %>%
    mutate(RR_diag = ifelse(N_diag == 0, NA, RR_diag))

  # Join to get RR(x,x) and RR(y,y) for each row
  df_rr <- df_rr %>%
    left_join(rr_diag, by = c("date", "x" = "age")) %>%
    rename(RR_xx = RR_diag) %>%
    left_join(rr_diag, by = c("date", "y" = "age")) %>%
    rename(RR_yy = RR_diag) %>%
    mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))

  return(df_rr)
}

df_age_time <- normalized_age_rr(df_age_time)

# Define age group order (from age_time_RR_analysis.R aggregation)
age_order <- c("0-5y", "6-17y", "18-25y", "26-45y", "46-65y", "65-80y", "81y+")

# Set factor levels for both plots
df_age_time <- df_age_time %>%
  mutate(
    x = factor(x, levels = age_order),
    y = factor(y, levels = age_order)
  )

# Plot 1: RR over time - include diagonal (all pairs)
plot_data_rr <- df_age_time

p1 <- ggplot(plot_data_rr, aes(x = date, y = RR, color = y, group = y)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", alpha = 0.5) +
  facet_wrap(~ x, scales = "free_y", nrow = 2, ncol = 4) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  scale_y_log10() +
  labs(
    title = "Relative Risk (RR) Over Time by Age Group",
    x = "Date",
    y = "RR (log scale)",
    color = "Paired Age Group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = c(0.875, 0.25),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.box = "vertical"
  )

ggsave(
  filename = paste0(fig_path, "age_RR_time_faceted.png"),
  plot = p1,
  width = 16,
  height = 8,
  dpi = 192
)

# Plot 2: nRR over time - exclude diagonal, place legend in 8th facet (bottom right)
plot_data_nrr <- df_age_time %>%
  filter(x != y)

p2 <- ggplot(plot_data_nrr, aes(x = date, y = nRR, color = y, group = y)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5) +
  facet_wrap(~ x, scales = "free_y", nrow = 2, ncol = 4) +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
  labs(
    title = "Normalized Relative Risk (nRR) Over Time by Age Group",
    x = "Date",
    y = "nRR",
    color = "Paired Age Group"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    legend.position = c(0.875, 0.25),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.box = "vertical"
  )

ggsave(
  filename = paste0(fig_path, "age_nRR_time_faceted.png"),
  plot = p2,
  width = 16,
  height = 8,
  dpi = 192
)

cat("Plots saved to:", fig_path, "\n")
