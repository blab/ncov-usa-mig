library(tidyverse)
library(data.table)
library(argparse)
library(igraph)
library(ggraph)
library(sf)

source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--sd_threshold', type = 'double', default = 3.0, help = 'Number of standard deviations for significance threshold')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
sd_threshold <- args$sd_threshold

select <- dplyr::select

# Read region data for coloring nodes by region
REGION_DATA <- fread("data/us_states_regions.csv")

# Load CAM shapefile and extract centroids
cam_map <- st_read("data/shp-files/cam-shp.gpkg", quiet = TRUE)
st_crs(cam_map) <- 3857  # Fix CRS
cam_map <- st_transform(cam_map, 4326)  # Transform to WGS84

# Calculate centroids for each state
state_centroids <- cam_map %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[, 1],
    lat = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  as_tibble() %>%
  select(state = NAME_En, lon, lat)

# Read the significant connections
fn_in <- paste0("results/", scenario, "/time_state/df_significant_connections_", sd_threshold, "sd.tsv")
sig_connections <- fread(fn_in)

# Create output directory for network plots
fn_out_path <- paste0("figs/", scenario, "/time_state_networks/")
dir.create(file.path(fn_out_path), showWarnings = FALSE, recursive = TRUE)

# Get unique time points
time_points <- unique(sig_connections$date) %>% sort()

cat("Creating network visualizations for", length(time_points), "time points...\n")

# Create a network plot for each time point
for(i in 1:length(time_points)){
  current_date <- time_points[i]

  # Filter to current time point
  current_edges <- sig_connections %>%
    filter(date == current_date) %>%
    select(x, y, nRR, z_score_x, z_score_y, N_pairs)

  # Skip if no edges
  if(nrow(current_edges) == 0){
    next
  }

  # Create graph from edge list
  g <- graph_from_data_frame(current_edges, directed = FALSE)

  # Add node attributes (region, country)
  node_names <- V(g)$name
  node_data <- REGION_DATA %>%
    filter(state %in% node_names) %>%
    select(state, bea_reg, country)

  # Match node order
  node_data <- node_data[match(node_names, node_data$state), ]

  V(g)$region <- node_data$bea_reg
  V(g)$country <- node_data$country

  # Add geographic coordinates for layout
  node_coords <- state_centroids %>%
    filter(state %in% node_names)

  # Match node order and create layout matrix
  node_coords <- node_coords[match(node_names, node_coords$state), ]
  layout_matrix <- as.matrix(node_coords[, c("lon", "lat")])

  # Edge weights based on minimum z-score (most conservative)
  E(g)$weight <- pmin(abs(current_edges$z_score_x), abs(current_edges$z_score_y))
  E(g)$n_pairs <- current_edges$N_pairs

  # Create the plot using standardized color schemes and geographic layout
  p <- ggraph(g, layout = layout_matrix) +
    geom_edge_link(aes(width = weight, alpha = weight), color = "gray40") +
    geom_node_point(aes(color = region), size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    region_color_scale() +  # Use standardized region colors from color_schemes.R
    scale_edge_width_continuous(range = c(0.5, 3), name = "Min |Z-score|") +
    scale_edge_alpha_continuous(range = c(0.3, 0.9), guide = "none") +
    labs(title = paste0("Significant Inter-State Connections: ", current_date),
         subtitle = paste0("Threshold: ", sd_threshold, " SD from out-of-region mean (both directions)\n",
                          "Edges: ", nrow(current_edges), " | Nodes: ", vcount(g))) +
    theme_graph() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    )

  # Save plot
  fn_out <- paste0(fn_out_path, "network_", format(current_date, "%Y-%m-%d"), ".jpg")
  ggsave(fn_out, plot = p, width = 12, height = 8, dpi = 300,create.dir=TRUE)

  cat("  [", i, "/", length(time_points), "] Saved:", fn_out, "\n")
}

cat("\nDone! Network plots saved to:", fn_out_path, "\n")

# Create summary statistics
cat("\n=== Network Summary Statistics by Time Point ===\n")
network_stats <- sig_connections %>%
  group_by(date) %>%
  summarise(
    n_edges = n(),
    n_nodes = n_distinct(c(x, y)),
    mean_nRR = mean(nRR, na.rm = TRUE),
    median_nRR = median(nRR, na.rm = TRUE),
    mean_z_x = mean(abs(z_score_x), na.rm = TRUE),
    mean_z_y = mean(abs(z_score_y), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date)

print(network_stats)

# Save summary stats
write_tsv(network_stats,
          file = paste0(fn_out_path, "network_summary_stats.tsv"))

cat("\n=== Network Analysis Complete ===\n")
