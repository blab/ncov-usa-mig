#!/usr/bin/env Rscript
# Correlation between age-specific RR and empirical contact matrices
# Compares RR matrices (within-country) against contact survey data
# Contact matrices from: https://github.com/mobs-lab/mixing-patterns
# Contact rates are converted to RR for direct comparison

library(tidyverse)
library(argparse)
library(duckdb)
library(dbplyr)

source("scripts/color_schemes.R")
source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  parser$add_argument('--exclude_duplicates', type = 'logical', default = FALSE,
                      help = "Whether to exclude possible duplicate pairs")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
exclude_duplicates <- args$exclude_duplicates

AGE_CAP <- 80  # Cap analyses at age 80

# --- Load contact matrices and convert to RR ---
# Contact RR: RR(i,j) = (C_ij * N_total) / (C_i. * C_.j)
# where C_i. = sum of row i (total contacts for age i)
load_contact_as_rr <- function(fn) {
  mat <- as.matrix(read.csv(fn, header = FALSE))
  mat <- mat[1:AGE_CAP, 1:AGE_CAP]

  # Row and column marginals
  row_sums <- rowSums(mat)
  col_sums <- colSums(mat)
  total <- sum(mat)

  # Compute RR matrix
  rr_mat <- matrix(0, nrow = AGE_CAP, ncol = AGE_CAP)
  for(i in 1:AGE_CAP) {
    for(j in 1:AGE_CAP) {
      rr_mat[i, j] <- (mat[i, j] * total) / (row_sums[i] * col_sums[j])
    }
  }

  # Convert to long format
  df <- expand.grid(
    x = sprintf("%02dy", 0:(AGE_CAP - 1)),
    y = sprintf("%02dy", 0:(AGE_CAP - 1)),
    stringsAsFactors = FALSE
  )
  df$contact_RR <- as.vector(rr_mat)
  return(df)
}

df_contact_us <- load_contact_as_rr("data/us_contact_matrix.csv")
df_contact_ca <- load_contact_as_rr("data/canada_contact_matrix.csv")

message("Contact matrices loaded and converted to RR.")

# --- Calculate country-specific single-year age RR ---
fn_db <- paste0("db_files/db_", scenario, ".duckdb")
con <- DBI::dbConnect(duckdb(), fn_db)

# Create country pairs binding
pairs_country <- con %>%
  bind_pairs_exp("country", exclude_duplicates = exclude_duplicates) %>%
  rename(country.x = x, country.y = y)

# Calculate RR for each country (within-country pairs only)
calculate_country_rr <- function(con, country_name) {
  message(paste0("Calculating RR for ", country_name, "..."))

  pairs_age <- con %>%
    bind_pairs_exp("age_class", exclude_duplicates = exclude_duplicates) %>%
    left_join(pairs_country, join_by(strain_1, strain_2)) %>%
    filter(country.x == country_name & country.y == country_name) %>%
    collect()

  rr <- pairs_age %>%
    calculate_rr_matrix()

  return(rr)
}

df_rr_us <- calculate_country_rr(con, "USA")
df_rr_ca <- calculate_country_rr(con, "Canada")

DBI::dbDisconnect(con, shutdown = TRUE)
message("RR calculations complete.")

# Save country-specific RR
fn_rr_us <- paste0("results/", scenario, "/df_RR_by_age_USA.tsv")
fn_rr_ca <- paste0("results/", scenario, "/df_RR_by_age_Canada.tsv")
write_tsv(df_rr_us, fn_rr_us)
write_tsv(df_rr_ca, fn_rr_ca)
message("Country RR matrices saved.")

# --- Cap RR at AGE_CAP ---
cap_rr <- function(df_rr) {
  df_rr %>%
    mutate(x_age = as.numeric(gsub("[^0-9]", "", x)),
           y_age = as.numeric(gsub("[^0-9]", "", y))) %>%
    filter(x_age < AGE_CAP & y_age < AGE_CAP) %>%
    select(-x_age, -y_age)
}

df_rr_us <- cap_rr(df_rr_us)
df_rr_ca <- cap_rr(df_rr_ca)

# --- Join RR with contact RR ---
df_corr_us <- df_rr_us %>%
  inner_join(df_contact_us, by = c("x", "y")) %>%
  mutate(country = "USA")

df_corr_ca <- df_rr_ca %>%
  inner_join(df_contact_ca, by = c("x", "y")) %>%
  mutate(country = "Canada")

df_corr <- bind_rows(df_corr_us, df_corr_ca) %>%
  filter(N_pairs > 0) %>%
  mutate(
    log_RR = log10(RR),
    log_contact_RR = log10(contact_RR)
  )

# Calculate Pearson correlations on log-transformed values
corr_stats <- df_corr %>%
  group_by(country) %>%
  summarize(
    r = cor(log_RR, log_contact_RR, method = "pearson", use = "complete.obs"),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r = %.2f", r))

message("\nCorrelation results:")
print(corr_stats)

# --- Plot ---
df_corr <- df_corr %>% mutate(same_age = (x == y))

p <- ggplot(df_corr, aes(x = contact_RR, y = RR, color = same_age)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(data = df_corr %>% filter(!same_age), alpha = 0.3, size = 0.5) +
  geom_point(data = df_corr %>% filter(same_age), alpha = 0.5, size = 1) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "firebrick1"),
                     labels = c("Different", "Same"), name = "Age") +
  geom_text(data = corr_stats,
            aes(x = 0.2, y = max(df_corr$RR, na.rm = TRUE), label = label),
            hjust = 0, vjust = 1, size = 2.5, inherit.aes = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  facet_wrap(~ country) +
  labs(
    x = expression(RR[Contact~Matrix]),
    y = expression(RR[Sequence])
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90"),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  )

fn_plot <- paste0("figs/", scenario, "/age_contact_correlation.png")
dir.create(dirname(fn_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(fn_plot, p, width = 7, height = 3, dpi = 300)
ggsave(sub("\\.png$", ".svg", fn_plot), p, width = 7, height = 3)
ggsave(sub("\\.png$", "_compact.svg", fn_plot), p, width = 5, height = 3)
message(paste0("Plot saved to ", fn_plot))

cat("\nAge-contact correlation analysis complete.\n")
