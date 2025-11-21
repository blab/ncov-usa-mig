# File: state_heatmap_clustered.R
# Author(s): Amin Bemanian
# Date: 2025-11-02
# Description: Creates a combined figure with hierarchical clustering dendrogram
#              on the left (y-axis) and state RR heatmap on the right, with
#              states ordered by clustering and colored by region
# Arguments:
#   --scenario: Scenario corresponding to data files (e.g. CAM_1000)
#   --dist: Cluster Distance Function Reciprocal (R) or Neg-Exponential (E), default = "E"

library(argparse)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggdendro)
library(cluster)
library(patchwork)
source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on',default = "CAM_1000")
  parser$add_argument('--dist', type = 'character', help = "Cluster Distance Function Reciprocal (R) or Neg-Exponential (E)", default = "E")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
dist_type <- args$dist

# Set region coloring
region_col <- "bea_reg"
label_title <- "BEA Regions"
get_region <- function(df) df[[region_col]]

# Load region data
df_regions <- read_csv("data/us_states_regions.csv") %>%
  mutate(hhs_reg = factor(hhs_reg,
                          levels = c("Boston (Region I)",
                                     "New York (Region II)",
                                     "Philadelphia (Region III)",
                                     "Atlanta (Region IV)",
                                     "Chicago (Region V)",
                                     "Dallas (Region VI)",
                                     "Kansas City (Region VII)",
                                     "Denver (Region VIII)",
                                     "San Francisco (Region IX)",
                                     "Seattle (Region X)",
                                     "Mexico",
                                     "Eastern Canada",
                                     "Western Canada"))) %>%
  mutate(census_div = factor(census_div,
                          levels = c("East North Central",
                                     "East South Central",
                                     "Middle Atlantic",
                                     "Mountain",
                                     "New England",
                                     "Pacific",
                                     "South Atlantic",
                                     "West North Central",
                                     "West South Central",
                                     "Mexico",
                                     "Eastern Canada",
                                     "Western Canada"))) %>%
  mutate(bea_reg = factor(bea_reg,
                             levels = c("Far West",
                                        "Great Lakes",
                                        "Mideast",
                                        "New England",
                                        "Plains",
                                        "Rocky Mountain",
                                        "Southeast",
                                        "Southwest",
                                        "Mexico",
                                        "Eastern Canada",
                                        "Western Canada"))) %>%
  mutate(mig_networks = factor(mig_networks))

# Load the full time period state RR matrix for clustering
fn_path_rr <- paste0("results/", scenario, "/time_state/")
state_rr_all <- fread(paste0(fn_path_rr, "df_state_rr_all.tsv"))

# Convert to matrix for clustering
state_rr_all_mat <- state_rr_all %>%
  select(c("x", "y", "RR")) %>%
  rename(name = x) %>%
  pivot_wider(names_from = y, values_from = RR) %>%
  arrange(name) %>%
  select(name, sort(setdiff(names(.), "name"))) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Create distance matrix
if (dist_type == 'E') {
  state_dist_mat <- exp(-state_rr_all_mat) %>% as.matrix() %>% as.dist()
  dist_full_name <- "Negative Exponential"
} else if (dist_type == 'R') {
  state_dist_mat <- 1 / (state_rr_all_mat) %>% as.matrix() %>% as.dist()
  dist_full_name <- "Reciprocal"
} else {
  stop("Unsupported dist_type: must be 'E' or 'R'")
}

# Perform hierarchical clustering (AGNES with ward method)
hc_agnes_full <- agnes(state_dist_mat, method = "ward", diss = TRUE) %>% as.hclust()

# Get dendrogram order
dendro_data_full <- dendro_data(hc_agnes_full, type = "rectangle")

# Extract state order from dendrogram (sort by x position which represents leaf order)
state_order <- dendro_data_full$labels %>%
  arrange(x) %>%
  pull(label)

# Create dendrogram plot for y-axis (vertical, on left side)
# Only show branches, no leaf points - labels will be on the heatmap
label_df <- dendro_data_full$labels %>%
  left_join(df_regions, join_by(label == state)) %>%
  mutate(region_label = get_region(.))

# Shift dendrogram segments to the left to avoid overlapping with labels
# The terminal branches end at y=0, shift them left so labels can be at x=0.5
# With scale_x_reverse, we need to ADD to push branches left
dendro_segments_shifted <- dendro_data_full$segments %>%
  mutate(y = y,
         yend = yend)

gg_dendro <- ggplot() +
  geom_segment(data = dendro_segments_shifted,
               aes(x = y, y = x, xend = yend, yend = xend)) +
  scale_x_reverse(expand = expansion(mult = c(0.15, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_void() +
  theme(
    plot.margin = unit(c(0.5, 0, 1, 0.2), "cm")
  )

# Load heatmap data
fn_rr <- paste("results/", scenario, "/df_RR_by_state.tsv", sep = "")
state_rr <- fread(fn_rr)

# Set RR bounds for heatmap
UB <- 2
LB <- 0.5


# Add region information for axis labeling
state_rr_labeled <- state_rr %>%
  left_join(df_regions %>% select(state, region = !!sym(region_col)),
            by = c("x" = "state")) %>%
  rename(x_region = region) %>%
  left_join(df_regions %>% select(state, region = !!sym(region_col)),
            by = c("y" = "state")) %>%
  rename(y_region = region)

# Create heatmap
AXIS_SIZE <- 7

# Create y-axis label data with region colors (in dendrogram order)
# Use 3 spaces as placeholder for uniform leaf indicators
y_label_data <- df_regions %>%
  filter(state %in% state_order) %>%
  mutate(state = factor(state, levels = state_order)) %>%
  arrange(state) %>%
  mutate(region = get_region(.),
         color = REGION_SCALE[as.character(region)],
         leaf_label = "               ",  # 3 spaces as placeholder
         y_pos = as.numeric(state),
         x_pos = 0.5)  # Position at left edge

# Reorder factors based on dendrogram (both axes use same order for symmetry)
state_rr_labeled$x <- factor(state_rr_labeled$x, levels = state_order)
state_rr_labeled$y <- factor(state_rr_labeled$y, levels = state_order)
state_rr_labeled <- state_rr_labeled %>% arrange(x)

state_heatmap <- state_rr_labeled %>%
  rowwise() %>%
  mutate(fill_RR = fill_bound(RR, LB, UB)) %>%
  ggplot(aes(x = y, y = x, fill = fill_RR)) +
  geom_tile() +
  # Add colored leaf indicators (3 spaces as placeholder)
  geom_label(data = y_label_data,
             aes(x = x_pos-0.2, y = state, label = leaf_label),
             fill = y_label_data$color,
             color = "white",
             size = 2.5,
             hjust = 1,
             linewidth = 0,
             inherit.aes = FALSE) +
  RR_log_grad(LB, UB) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0.5, 0.2, 1, 0.1), "cm")
  ) +
  labs(x = "State/Province",
       y = NULL,
       title = "Identical Sequence RR with Hierarchical Clustering") +
  coord_fixed(clip = "off")

# Combine dendrogram and heatmap with minimal spacing
combined_plot <- gg_dendro + state_heatmap +
  plot_layout(widths = c(1, 4))

# Save combined figure
dir.create(paste0("figs/", scenario), recursive = TRUE, showWarnings = FALSE)
fn_combined <- paste0("figs/", scenario, "/state_heatmap_clustered.jpg")

ggsave(fn_combined,
       plot = combined_plot,
       device = "jpeg",
       dpi = 192,
       width = 12,
       height = 10,
       create.dir = TRUE)

message(paste("Combined dendrogram + heatmap saved to:", fn_combined))
