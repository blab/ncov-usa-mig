#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)
library(lubridate)
library(lme4)
library(lmerTest)
library(splines)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# Read the time-stratified RR data
fn_rr_time <- paste0("results/", scenario, "/time_age/df_RR_by_school_state_time.tsv")
df_rr <- read_tsv(fn_rr_time, show_col_types = FALSE)

# Read the school share data
fn_school <- "data/state_school_share.csv"
df_school <- read_csv(fn_school, show_col_types = FALSE)

# Convert dates and add year-month column
df_school <- df_school %>%
  mutate(
    date = as.Date(date),
    year_month = floor_date(date, "month")
  )

# ------------------------------------------------------------------
# Trajectory Classification
# ------------------------------------------------------------------

# Define academic year period (September 2020 - May 2021)
academic_year_start <- as.Date("2020-09-01")
academic_year_end <- as.Date("2021-05-31")

# Filter to academic year period
df_school_ay <- df_school %>%
  filter(date >= academic_year_start & date <= academic_year_end)

# Calculate trajectory metrics for each state
state_trajectories <- df_school_ay %>%
  group_by(state) %>%
  arrange(date) %>%
  summarize(
    start_inperson = first(share_inperson),
    end_inperson = last(share_inperson),
    min_inperson = min(share_inperson, na.rm = TRUE),
    max_inperson = max(share_inperson, na.rm = TRUE),
    mean_inperson = mean(share_inperson, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(
    trajectory_type = case_when(
      start_inperson >= 0.8 & min_inperson >= 0.8 ~ "Persistent High",
      start_inperson <= 0.2 & max_inperson <= 0.2 ~ "Persistent Low",
      start_inperson <= 0.2 & end_inperson >= 0.8 ~ "Rising In-Person",
      TRUE ~ "Other"
    )
  )

# Print classification results
cat("\n=== STATE TRAJECTORY CLASSIFICATIONS ===\n\n")

for(traj_type in c("Persistent High", "Persistent Low", "Rising In-Person", "Other")) {
  states_in_category <- state_trajectories %>%
    filter(trajectory_type == traj_type) %>%
    arrange(state) %>%
    pull(state)

  cat(sprintf("%s (%d states):\n", traj_type, length(states_in_category)))
  if(length(states_in_category) > 0) {
    cat(paste("  -", states_in_category, collapse = "\n"), "\n")
  }
  cat("\n")
}

# Save trajectory classifications
fn_trajectories_class <- paste0("results/", scenario, "/time_age/state_trajectory_classifications.tsv")
write_tsv(state_trajectories, fn_trajectories_class)
message(paste0("Trajectory classifications saved to ", fn_trajectories_class))

# Join trajectory classification back to full data
df_school <- df_school %>%
  left_join(state_trajectories %>% select(state, trajectory_type), by = "state")

# Create faceted plot by trajectory type for the three main categories
df_school_main_trajectories <- df_school %>%
  filter(trajectory_type %in% c("Persistent High", "Persistent Low", "Rising In-Person")) %>%
  mutate(
    trajectory_type = factor(trajectory_type,
                            levels = c("Persistent High", "Rising In-Person", "Persistent Low"))
  )

# Define colors for trajectory types
trajectory_colors <- c(
  "Persistent High" = "#2E7D32",    # Dark green
  "Persistent Low" = "#C62828",     # Dark red
  "Rising In-Person" = "#1565C0"    # Dark blue
)

p_trajectories <- ggplot(df_school_main_trajectories,
                         aes(x = date, y = share_inperson,
                             group = state, color = trajectory_type)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +
  facet_wrap(~ trajectory_type, ncol = 1) +
  scale_color_manual(values = trajectory_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    x = "Date",
    y = "In-Person Instruction Share",
    title = "State In-Person Instruction Trajectories by Classification",
    subtitle = "Academic Year 2020-2021 (Sept 2020 - May 2021)",
    color = "Trajectory Type"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

fn_trajectories_plot <- paste0("figs/", scenario, "/age_time/school_trajectories_by_type.jpg")
ggsave(fn_trajectories_plot, p_trajectories, width = 10, height = 10, dpi = 300)
message(paste0("Trajectory plot saved to ", fn_trajectories_plot))

# Create heatmap version of trajectories
# Calculate the date interval for proper tile width (no gaps)
date_interval <- as.numeric(diff(sort(unique(df_school_main_trajectories$date)))[1])

p_trajectories_heatmap <- ggplot(df_school_main_trajectories,
                                  aes(x = date, y = state, fill = share_inperson)) +
  geom_tile(width = date_interval) +
  facet_grid(trajectory_type ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#C62828",      # Dark red (virtual)
    mid = "#FDD835",      # Yellow (hybrid)
    high = "#2E7D32",     # Dark green (in-person)
    midpoint = 0.5,
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent_format(accuracy = 1)
  ) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = "Date",
    y = "State",
    fill = "In-Person\nShare",
    title = "State In-Person Instruction Heatmap by Trajectory Classification",
    subtitle = "Academic Year 2020-2021 (Sept 2020 - May 2021)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

fn_trajectories_heatmap <- paste0("figs/", scenario, "/age_time/school_trajectories_heatmap.jpg")
ggsave(fn_trajectories_heatmap, p_trajectories_heatmap, width = 8, height = 6, dpi = 300)
message(paste0("Trajectory heatmap saved to ", fn_trajectories_heatmap))

# ------------------------------------------------------------------

# Filter RR data to within-group transmission for Primary and Secondary School
df_rr_school <- df_rr %>%
  filter(x == y) %>%
  filter(x %in% c("Primary School", "Secondary School")) %>%
  mutate(
    date = as.Date(date),
    year_month = floor_date(date, "month")
  ) %>%
  # Average RRs if there are multiple in the same state-month-age group
  group_by(state, x, year_month) %>%
  summarize(
    RR = mean(RR, na.rm = TRUE),
    nRR = mean(nRR, na.rm = TRUE),
    nRR_fixed = mean(nRR_fixed, na.rm = TRUE),
    N_pairs = sum(N_pairs),
    .groups = "drop"
  )

# Join RR data with school share data
df_combined <- df_rr_school %>%
  inner_join(df_school, by = c("state", "year_month" = "year_month"))

# ------------------------------------------------------------------
# 1. Time series plot of school shares
# ------------------------------------------------------------------

# Pivot longer for faceting
df_school_long <- df_school %>%
  pivot_longer(
    cols = c(share_inperson, share_hybrid, share_virtual),
    names_to = "share_type",
    values_to = "share_value"
  ) %>%
  mutate(
    share_type = factor(share_type,
                       levels = c("share_inperson", "share_hybrid", "share_virtual"),
                       labels = c("In-Person", "Hybrid", "Virtual"))
  )

p_timeseries <- ggplot(df_school_long, aes(x = date, y = share_value, group = state)) +
  geom_line(alpha = 0.3, color = "gray50") +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               color = "red", linewidth = 1.5) +
  facet_wrap(~ share_type, ncol = 1, scales = "free_x") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    x = "Date",
    y = "School Share",
    title = "School Instruction Modality Over Time by State",
    subtitle = "Gray lines = individual states, Red line = mean across states"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11)
  )

fn_timeseries <- paste0("figs/", scenario, "/age_time/school_share_timeseries.jpg")
dir.create(dirname(fn_timeseries), recursive = TRUE, showWarnings = FALSE)
ggsave(fn_timeseries, p_timeseries, width = 10, height = 10, dpi = 300)
message(paste0("Time series plot saved to ", fn_timeseries))

# ------------------------------------------------------------------
# 2. Correlation plots (2x3 grid)
# ------------------------------------------------------------------

# Prepare data for correlation plots
# Apply floor at 0.1 * minimum RR
rr_floor <- 0.1 * min(df_combined$RR, na.rm = TRUE)
df_corr <- df_combined %>%
  mutate(
    RR_floored = pmax(RR, rr_floor),
    log_RR = log(RR_floored)
  ) %>%
  pivot_longer(
    cols = c(share_inperson, share_hybrid, share_virtual),
    names_to = "share_type",
    values_to = "share_value"
  ) %>%
  mutate(
    share_type = factor(share_type,
                       levels = c("share_inperson", "share_hybrid", "share_virtual"),
                       labels = c("In-Person", "Hybrid", "Virtual")),
    x = factor(x, levels = c("Primary School", "Secondary School"))
  )

# Calculate Spearman correlations and p-values for annotations
df_corr_stats <- df_corr %>%
  group_by(x, share_type) %>%
  summarize(
    cor_test = list(cor.test(log_RR, share_value, method = "spearman")),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    rho = map_dbl(cor_test, ~ .x$estimate),
    p_value = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  mutate(
    label = sprintf("Spearman rho = %.2f", rho, p_value, n)
  ) %>%
  select(-cor_test)

p_corr <- ggplot(df_corr, aes(x = share_value, y = log_RR)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, color = "firebrick", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = df_corr_stats,
            aes(x = 0.05, y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 3.5) +
  facet_grid(share_type ~ x, scales = "free") +
  labs(
    x = "School Share",
    y = "log(Within-Group Relative Risk)",
    title = "Correlation Between School Instruction Modality and Within-Group Transmission"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10)
  )

fn_corr <- paste0("figs/", scenario, "/age_time/school_share_RR_correlation.jpg")
ggsave(fn_corr, p_corr, width = 10, height = 10, dpi = 300)
message(paste0("Correlation plot saved to ", fn_corr))

# ------------------------------------------------------------------
# 3. Correlation plots for nRR_fixed (adult baseline normalized)
# ------------------------------------------------------------------

df_corr_nrr_fixed <- df_combined %>%
  pivot_longer(
    cols = c(share_inperson, share_hybrid, share_virtual),
    names_to = "share_type",
    values_to = "share_value"
  ) %>%
  mutate(
    share_type = factor(share_type,
                       levels = c("share_inperson", "share_hybrid", "share_virtual"),
                       labels = c("In-Person", "Hybrid", "Virtual")),
    x = factor(x, levels = c("Primary School", "Secondary School"))
  )

df_corr_stats_nrr_fixed <- df_corr_nrr_fixed %>%
  group_by(x, share_type) %>%
  summarize(
    cor_test = list(cor.test(nRR_fixed, share_value, method = "spearman")),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    rho = map_dbl(cor_test, ~ .x$estimate),
    p_value = map_dbl(cor_test, ~ .x$p.value)
  ) %>%
  mutate(
    label = sprintf("Spearman rho = %.2f", rho, p_value, n)
  ) %>%
  select(-cor_test)

p_corr_nrr_fixed <- ggplot(df_corr_nrr_fixed, aes(x = share_value, y = nRR_fixed)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = df_corr_stats_nrr_fixed,
            aes(x = 0.05, y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 3.5) +
  facet_grid(share_type ~ x, scales = "free") +
  coord_cartesian(ylim=c(-0.5,100)) +
  labs(
    x = "School Share",
    y = "Within-Group nRR",
    title = "Correlation Between School Modality and Adult-Normalized Within-Group Transmission"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10)
  )

fn_corr_nrr_fixed <- paste0("figs/", scenario, "/age_time/school_share_nRR_fixed_correlation.jpg")
ggsave(fn_corr_nrr_fixed, p_corr_nrr_fixed, width = 10, height = 10, dpi = 300)
message(paste0("nRR_fixed correlation plot saved to ", fn_corr_nrr_fixed))

# Save combined data for further analysis
fn_combined <- paste0("results/", scenario, "/time_age/df_RR_school_share_combined.tsv")
write_tsv(df_combined, fn_combined)
message(paste0("Combined data saved to ", fn_combined))

nrr_share_model <- lmer(nRR_fixed  ~  share_hybrid + bs(date)  + (1|state),
                        data=df_combined %>% filter(x == "Secondary School"))
summary(nrr_share_model)
ranova(nrr_share_model)
AIC(nrr_share_model)
