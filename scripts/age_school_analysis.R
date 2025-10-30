#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)
library(lme4)
library(car)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# Read the academic year RR data
fn_school_state_ay <- paste0("results/", scenario, "/time_age/df_RR_by_school_state_ay.tsv")
df_ay <- read_tsv(fn_school_state_ay, show_col_types = FALSE)

# Filter to only include school categories for X (facet), but keep all y categories
df_plot <- df_ay %>%
  filter(x %in% c("Pre-School", "Primary School", "Secondary School")) %>%
  # Set factor levels for ordering
  mutate(
    x = factor(x, levels = c("Pre-School", "Primary School", "Secondary School")),
    y = factor(y, levels = c("Pre-School", "Primary School", "Secondary School", "Adult", "Seniors", "NA")),
    academic_year = factor(academic_year,
                          levels = c("2020-2021", "2021-2022", "2022-2023", "2023-2024"))
  )

# Create the faceted boxplot
p <- ggplot(df_plot, aes(x = academic_year, y = RR, fill = y)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~ x, ncol = 3, scales = "free_x") +
  scale_y_log10(
    limits = c(0.01, 100),
    breaks = c(0.01, 0.1, 1, 10, 100),
    labels = c("0.01", "0.1", "1", "10", "100")
  ) +
  labs(
    x = "Academic Year",
    y = "Relative Risk (RR)",
    fill = "Age Group",
    title = "Distribution of RR Across States by School Stage and Academic Year"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

# Save the plot
fn_out <- paste0("figs/", scenario, "/age_time/age_school_ay_boxplot.jpg")
dir.create(dirname(fn_out), recursive = TRUE, showWarnings = FALSE)
ggsave(fn_out, p, width = 12, height = 6, dpi = 300)

message(paste0("Plot saved to ", fn_out))

# ------------------------------------------------------------------
# Additional plots for nRR and nRR_fixed
# ------------------------------------------------------------------

# Plot for nRR (diagonal normalized)
p_nrr <- ggplot(df_plot, aes(x = academic_year, y = nRR, fill = y)) +
  geom_boxplot() +
  facet_wrap(~ x, ncol = 3, scales = "free_x") +
  labs(
    x = "Academic Year",
    y = "Normalized Relative Risk (nRR)",
    fill = "Age Group",
    title = "Distribution of nRR (Diagonal Normalized) Across States by School Stage and Academic Year"
  ) +
  coord_cartesian(ylim=c(-0.5,5)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )

fn_nrr <- paste0("figs/", scenario, "/age_time/age_school_ay_nRR_boxplot.jpg")
ggsave(fn_nrr, p_nrr, width = 12, height = 6, dpi = 300)
message(paste0("nRR plot saved to ", fn_nrr))

# Plot for nRR_fixed (adult baseline normalized)
p_nrr_fixed <- ggplot(df_plot, aes(x = academic_year, y = nRR_fixed, fill = y)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  facet_wrap(~ x, ncol = 3, scales = "free_x") +
  labs(
    x = "Academic Year",
    y = "Fixed Baseline Normalized RR (nRR_fixed)",
    fill = "Age Group",
    title = "Distribution of nRR_fixed (Adult Baseline) Across States by School Stage and Academic Year"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim=c(-0.5,20))

fn_nrr_fixed <- paste0("figs/", scenario, "/age_time/age_school_ay_nRR_fixed_boxplot.jpg")
ggsave(fn_nrr_fixed, p_nrr_fixed, width = 12, height = 6, dpi = 300)
message(paste0("nRR_fixed plot saved to ", fn_nrr_fixed))

# Repeated measures ANOVA analysis
# Focus on within-group transmission (x == y) to test temporal trends
# State is the repeated measures factor

# Prepare data for ANOVA - filter to within-group RR (x == y) only
df_anova <- df_ay %>%
  filter(x %in% c("Pre-School", "Primary School", "Secondary School")) %>%
  filter(x == y) %>%  # Within-group transmission only
  mutate(log_RR = log(RR)) %>%
  filter(is.finite(log_RR)) %>%
  mutate(
    x = factor(x, levels = c("Pre-School", "Primary School", "Secondary School")),
    academic_year = factor(academic_year, levels = c("2020-2021", "2021-2022", "2022-2023", "2023-2024")),
    state = factor(state)
  )

# Fit separate models for each child age group to test temporal trends
child_groups <- c("Pre-School", "Primary School", "Secondary School")
anova_results <- list()

for (child_group in child_groups) {
  message(paste0("\n=== ANOVA for ", child_group, " (within-group transmission) ==="))

  # Filter to current child group
  df_subset <- df_anova %>% filter(x == child_group)

  # Skip if insufficient data
  if (nrow(df_subset) < 10) {
    message(paste0("Insufficient data for ", child_group))
    next
  }

  # Fit mixed effects model: log(RR) ~ academic_year + (1|state)
  # Tests whether within-group transmission risk changes over academic years
  # (1|state) = random intercept for state (repeated measures)
  model <- lmer(log_RR ~ academic_year + (1|state), data = df_subset)

  # Get Type II ANOVA table
  anova_table <- Anova(model, type = "II")

  # Store results
  anova_results[[child_group]] <- list(
    model = model,
    anova = anova_table,
    n_obs = nrow(df_subset),
    n_states = length(unique(df_subset$state)),
    summary = summary(model)
  )

  # Print results
  print(summary(model))
  print(anova_table)
  message(paste0("N observations: ", nrow(df_subset)))
  message(paste0("N states: ", length(unique(df_subset$state))))
}

# Also fit combined model to test for interaction between age group and time
message("\n=== Combined model: testing age group × year interaction ===")
model_combined <- lmer(log_RR ~ x * academic_year + (1|state), data = df_anova)
anova_combined <- Anova(model_combined, type = "II")
print(summary(model_combined))
print(anova_combined)

anova_results[["combined"]] <- list(
  model = model_combined,
  anova = anova_combined,
  n_obs = nrow(df_anova),
  n_states = length(unique(df_anova$state)),
  summary = summary(model_combined)
)

# Save ANOVA results
fn_anova_out <- paste0("results/", scenario, "/time_age/anova_school_ay_within_group_results.rds")
saveRDS(anova_results, fn_anova_out)
message(paste0("\nANOVA results saved to ", fn_anova_out))


# ------------------------------------------------------------------
# ANOVA for nRR_fixed (adult baseline normalized)
# ------------------------------------------------------------------

message("\n\n========================================")
message("ANOVA RESULTS FOR nRR_fixed (ADULT BASELINE NORMALIZED)")
message("========================================")

df_anova_nrr_fixed <- df_ay %>%
  filter(x %in% c("Pre-School", "Primary School", "Secondary School")) %>%
  filter(x == y) %>%
  mutate(log_nRR_fixed = log(nRR_fixed)) %>%
  filter(is.finite(log_nRR_fixed)) %>%
  mutate(
    x = factor(x, levels = c("Pre-School", "Primary School", "Secondary School")),
    academic_year = factor(academic_year, levels = c("2020-2021", "2021-2022", "2022-2023", "2023-2024")),
    state = factor(state)
  )

anova_results_nrr_fixed <- list()

for (child_group in child_groups) {
  message(paste0("\n=== ANOVA for ", child_group, " (within-group nRR_fixed) ==="))

  df_subset <- df_anova_nrr_fixed %>% filter(x == child_group)

  if (nrow(df_subset) < 10) {
    message(paste0("Insufficient data for ", child_group))
    next
  }

  model <- lmer(log_nRR_fixed ~ academic_year + (1|state), data = df_subset)
  anova_table <- Anova(model, type = "II")

  anova_results_nrr_fixed[[child_group]] <- list(
    model = model,
    anova = anova_table,
    n_obs = nrow(df_subset),
    n_states = length(unique(df_subset$state)),
    summary = summary(model)
  )

  print(summary(model))
  print(anova_table)
  message(paste0("N observations: ", nrow(df_subset)))
  message(paste0("N states: ", length(unique(df_subset$state))))
}

message("\n=== Combined model: testing age group × year interaction (nRR_fixed) ===")
model_combined_nrr_fixed <- lmer(log_nRR_fixed ~ x * academic_year + (1|state), data = df_anova_nrr_fixed)
anova_combined_nrr_fixed <- Anova(model_combined_nrr_fixed, type = "II")
print(summary(model_combined_nrr_fixed))
print(anova_combined_nrr_fixed)

anova_results_nrr_fixed[["combined"]] <- list(
  model = model_combined_nrr_fixed,
  anova = anova_combined_nrr_fixed,
  n_obs = nrow(df_anova_nrr_fixed),
  n_states = length(unique(df_anova_nrr_fixed$state)),
  summary = summary(model_combined_nrr_fixed)
)

fn_anova_nrr_fixed <- paste0("results/", scenario, "/time_age/anova_school_ay_nRR_fixed_results.rds")
saveRDS(anova_results_nrr_fixed, fn_anova_nrr_fixed)
message(paste0("\nnRR_fixed ANOVA results saved to ", fn_anova_nrr_fixed))


