# File: age_cluster.R
# Author: Amin Bemanian
# Date: 01/26/26
# Description: Hierarchical clustering analysis of age groups based on RR transmission patterns

library(tidyverse)
library(argparse)
library(ggplot2)
library(ggdendro)
library(patchwork)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', default = "CAM_1000", type = 'character',
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# Read the RR data
df_rr <- read_tsv(paste0("results/", scenario, "/df_RR_by_age_class.tsv"),
                  show_col_types = FALSE)

# Create RR matrix
rr_matrix <- df_rr %>%
  select(x, y, RR) %>%
  pivot_wider(names_from = y, values_from = RR) %>%
  column_to_rownames("x") %>%
  as.matrix()

cat("RR matrix dimensions:", dim(rr_matrix), "\n")
cat("Range of RR values:", round(range(rr_matrix, na.rm=TRUE), 2), "\n\n")

# ============================================================================
# Main clustering (5 clusters)
# ============================================================================

# Use Euclidean distance on the RR profiles
dist_matrix <- dist(rr_matrix, method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")

# Cut tree into 5 clusters
n_clusters <- 5
clusters_5 <- cutree(hc, k = n_clusters)

cat("=== MAIN CLUSTERING (5 CLUSTERS) ===\n\n")

for(i in 1:n_clusters) {
  ages_in_cluster <- names(clusters_5[clusters_5 == i])

  within_rr <- df_rr %>%
    filter(x %in% ages_in_cluster, y %in% ages_in_cluster) %>%
    summarize(mean_RR = mean(RR, na.rm = TRUE)) %>%
    pull(mean_RR)

  between_rr <- df_rr %>%
    filter(x %in% ages_in_cluster, !(y %in% ages_in_cluster)) %>%
    summarize(mean_RR = mean(RR, na.rm = TRUE)) %>%
    pull(mean_RR)

  cat("Cluster", i, "(n=", length(ages_in_cluster), "):\n")
  cat("  Ages:", paste(ages_in_cluster, collapse = ", "), "\n")
  cat("  Within-cluster RR:", round(within_rr, 2),
      "| Between-cluster RR:", round(between_rr, 2),
      "| Ratio:", round(within_rr/between_rr, 2), "\n\n")
}

# ============================================================================
# Sub-clustering within large cluster
# ============================================================================

cluster_sizes <- table(clusters_5)
large_cluster_id <- which.max(cluster_sizes)
large_cluster_ages <- names(clusters_5[clusters_5 == large_cluster_id])

cat("=== SUB-CLUSTERING (Working-age group:",
    min(large_cluster_ages), "-", max(large_cluster_ages), ") ===\n\n")

# Extract submatrix for the large cluster
submatrix <- rr_matrix[large_cluster_ages, large_cluster_ages]

# Perform sub-clustering
subdist <- dist(submatrix, method = "euclidean")
subhc <- hclust(subdist, method = "ward.D2")
subclusters <- cutree(subhc, k = 6)

for(i in 1:6) {
  ages_in_subcluster <- names(subclusters[subclusters == i])

  within_rr <- df_rr %>%
    filter(x %in% ages_in_subcluster, y %in% ages_in_subcluster) %>%
    summarize(mean_RR = mean(RR, na.rm = TRUE)) %>%
    pull(mean_RR)

  other_subcluster_ages <- large_cluster_ages[!large_cluster_ages %in% ages_in_subcluster]
  between_rr <- df_rr %>%
    filter(x %in% ages_in_subcluster, y %in% other_subcluster_ages) %>%
    summarize(mean_RR = mean(RR, na.rm = TRUE)) %>%
    pull(mean_RR)

  cat("Sub-cluster", i, "(n=", length(ages_in_subcluster), "):\n")
  cat("  Age range:", min(ages_in_subcluster), "-", max(ages_in_subcluster), "\n")
  cat("  Within-subcluster RR:", round(within_rr, 2),
      "| Ratio:", round(within_rr/between_rr, 2), "\n\n")
}

# ============================================================================
# Visualizations
# ============================================================================

fig_path <- paste0("figs/", scenario, "/age_cluster/")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

# Color palette for clusters
cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

# 1. Main dendrogram with cluster colors
dend_data <- dendro_data(hc, type = "rectangle")
dend_segments <- dend_data$segments

# Add cluster assignments to labels
dend_labels <- dend_data$labels %>%
  mutate(cluster = clusters_5[as.character(label)],
         cluster_color = cluster_colors[cluster])

p1 <- ggplot() +
  geom_segment(data = dend_segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dend_labels,
            aes(x = x, y = y, label = label, color = cluster_color),
            hjust = 1, angle = 90, size = 2) +
  scale_color_identity() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Hierarchical Clustering of Age Groups",
       subtitle = "5 clusters based on RR transmission patterns (Ward's method)",
       y = "Height (Euclidean distance)", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"))

ggsave(paste0(fig_path, "main_dendrogram.jpg"), p1,
       width = 16, height = 8, dpi = 300)

# 2. Sub-clustering dendrogram
subdend_data <- dendro_data(subhc, type = "rectangle")
subdend_segments <- subdend_data$segments

subcluster_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#35FAD9")

subdend_labels <- subdend_data$labels %>%
  mutate(subcluster = subclusters[as.character(label)],
         subcluster_color = subcluster_colors[subcluster])

p2 <- ggplot() +
  geom_segment(data = subdend_segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = subdend_labels,
            aes(x = x, y = y, label = label, color = subcluster_color),
            hjust = 1, angle = 90, size = 2.5) +
  scale_color_identity() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Sub-clustering: Pediatric and Working-Age Population (1-68y)",
       subtitle = "6 sub-clusters within the largest cluster",
       y = "Height (Euclidean distance)", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"))

ggsave(paste0(fig_path, "subcluster_dendrogram.jpg"), p2,
       width = 16, height = 8, dpi = 300)

# 3. Cluster summary barplot
cluster_summary <- tibble(
  cluster = 1:n_clusters,
  ages = sapply(1:n_clusters, function(i) {
    ages <- names(clusters_5[clusters_5 == i])
    paste0(min(ages), "-", max(ages))
  }),
  n = as.numeric(table(clusters_5)),
  within_rr = sapply(1:n_clusters, function(i) {
    ages <- names(clusters_5[clusters_5 == i])
    df_rr %>%
      filter(x %in% ages, y %in% ages) %>%
      summarize(mean_RR = mean(RR, na.rm = TRUE)) %>%
      pull(mean_RR)
  }),
  between_rr = sapply(1:n_clusters, function(i) {
    ages <- names(clusters_5[clusters_5 == i])
    df_rr %>%
      filter(x %in% ages, !(y %in% ages)) %>%
      summarize(mean_RR = mean(RR, na.rm = TRUE)) %>%
      pull(mean_RR)
  })
) %>%
  mutate(ratio = within_rr / between_rr,
         label = paste0("Cluster ", cluster, "\n(", ages, ")\nn=", n))

p3 <- ggplot(cluster_summary, aes(x = reorder(label, -cluster), y = ratio, fill = factor(cluster))) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Within-cluster / Between-cluster RR Ratio",
       subtitle = "Higher ratios indicate stronger preferential mixing within age group",
       x = "", y = "Within/Between RR Ratio") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold"))

ggsave(paste0(fig_path, "cluster_ratios.jpg"), p3,
       width = 10, height = 6, dpi = 300)


cat("\n=== Plots saved to:", fig_path, "===\n")
cat("  - main_dendrogram.jpg\n")
cat("  - subcluster_dendrogram.jpg\n")
cat("  - cluster_ratios.jpg\n")