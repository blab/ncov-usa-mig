library(tidyverse)
library(data.table)
library(argparse)
library(igraph)
library(ggraph)
library(sf)

source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--sd_threshold', type = 'double', default = 2, help = 'Number of standard deviations for significance threshold')
  parser$add_argument('--min_occurrences', type = 'integer', default = 3, help = 'Minimum number of time points a connection must appear in')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
sd_threshold <- args$sd_threshold
min_occurrences <- args$min_occurrences

select <- dplyr::select

# Read region data
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

# Create edge pairs (unordered, so A-B is same as B-A)
sig_connections <- sig_connections %>%
  mutate(
    edge_pair = paste(pmin(x, y), pmax(x, y), sep = "_")
  )

# Count how many times each edge appears across all time points
edge_counts <- sig_connections %>%
  group_by(edge_pair, x, y, region_x, region_y, country_x, country_y) %>%
  summarise(
    n_occurrences = n(),
    first_date = min(date),
    last_date = max(date),
    mean_nRR = mean(nRR, na.rm = TRUE),
    median_nRR = median(nRR, na.rm = TRUE),
    mean_z_x = mean(abs(z_score_x), na.rm = TRUE),
    mean_z_y = mean(abs(z_score_y), na.rm = TRUE),
    max_z = max(pmin(abs(z_score_x), abs(z_score_y))),
    .groups = "drop"
  ) %>%
  filter(n_occurrences >= min_occurrences) %>%
  arrange(desc(n_occurrences), desc(max_z))

cat("\n=== Conserved Connection Analysis ===\n")
cat("Connections appearing in", min_occurrences, "or more time points:", nrow(edge_counts), "\n")
cat("\nTop 10 most conserved connections:\n")
print(head(edge_counts %>% select(x, y, n_occurrences, first_date, last_date, mean_nRR), 10))

# Create edges for the conserved network
conserved_edges <- edge_counts %>%
  select(x, y, n_occurrences, mean_nRR, mean_z_x, mean_z_y, max_z)

# Create graph from edge list
g <- graph_from_data_frame(conserved_edges, directed = FALSE)

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

# Match node order
node_coords <- node_coords[match(node_names, node_coords$state), ]

# Add small jitter to prevent overlapping nodes in dense areas
# Use smaller jitter (0.3 degrees ~ 20 miles) to maintain geographic accuracy
set.seed(42)  # For reproducibility
jitter_amount <- 0.3
node_coords <- node_coords %>%
  mutate(
    lon = lon + runif(n(), -jitter_amount, jitter_amount),
    lat = lat + runif(n(), -jitter_amount, jitter_amount)
  )

# Create layout matrix
layout_matrix <- as.matrix(node_coords[, c("lon", "lat")])

# Edge weights based on number of occurrences
E(g)$weight <- conserved_edges$n_occurrences
E(g)$max_z <- conserved_edges$max_z

# Create the plot
p <- ggraph(g, layout = layout_matrix) +
  geom_edge_link(aes(width = weight, alpha = weight), color = "gray40") +
  geom_node_point(aes(color = region), size = 6) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3.5, fontface = "bold") +
  region_color_scale() +
  scale_edge_width_continuous(range = c(0.5, 4),
                               name = "# Time Points",
                               breaks = c(2, 5, 10, 15)) +
  scale_edge_alpha_continuous(range = c(0.4, 0.95), guide = "none") +
  labs(title = "Conserved Significant Inter-State Connections",
       subtitle = paste0("Connections appearing in ", min_occurrences, "+ time points (out of 19 quarters)\n",
                        "Total conserved connections: ", nrow(conserved_edges), " | States involved: ", vcount(g))) +
  theme_graph() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "right"
  )

# Save plot
fn_out_path <- paste0("figs/", scenario, "/time_state_networks/")
fn_out <- paste0(fn_out_path, "network_conserved_", min_occurrences, "plus.jpg")
ggsave(fn_out, plot = p, width = 14, height = 9, dpi = 300,create.dir=TRUE)

cat("\nConserved network plot saved to:", fn_out, "\n")

# Save the conserved edges data
fn_out_data <- paste0("results/", scenario, "/time_state/df_conserved_connections_", min_occurrences, "plus.tsv")
write_tsv(edge_counts, fn_out_data)

cat("Conserved connections data saved to:", fn_out_data, "\n")

# ============================================================================
# GEOGRAPHIC MAP VISUALIZATION
# ============================================================================
cat("\n=== Creating Geographic Map Overlay ===\n")

# Get unique states involved in conserved connections
states_involved <- unique(c(edge_counts$x, edge_counts$y))

# Prepare edge data for plotting
# Keep only unique edges (remove duplicates where x-y appears twice as x→y and y→x)
# Strategy: keep the edge where the western state (lower longitude) is the source
edges_for_plot <- edge_counts %>%
  left_join(state_centroids, by = c("x" = "state")) %>%
  rename(lon_x = lon, lat_x = lat) %>%
  left_join(state_centroids, by = c("y" = "state")) %>%
  rename(lon_y = lon, lat_y = lat) %>%
  # For each edge pair, keep only the one where lon_x <= lon_y (westernmost as source)
  # This ensures we only draw one arc per connection
  group_by(edge_pair) %>%
  slice_min(lon_x, n = 1, with_ties = FALSE) %>%
  ungroup()

# Generate simple geodesic-like arcs by adding curvature perpendicular to the line
# This approximates great circle behavior for moderate distances
generate_curved_path <- function(lon1, lat1, lon2, lat2, n = 50, curvature = 0.15) {
  # Parameterize the path from 0 to 1
  t <- seq(0, 1, length.out = n)

  # Linear interpolation
  lon_linear <- lon1 + t * (lon2 - lon1)
  lat_linear <- lat1 + t * (lat2 - lat1)

  # Add perpendicular offset (creates the arc)
  # Maximum offset at midpoint (t=0.5)
  offset <- curvature * sin(t * pi)

  # Perpendicular direction (rotate 90 degrees)
  dx <- lon2 - lon1
  dy <- lat2 - lat1
  len <- sqrt(dx^2 + dy^2)

  # Normalized perpendicular vector
  perp_x <- -dy / len
  perp_y <- dx / len

  # Apply offset
  lon_curved <- lon_linear + offset * perp_x * len
  lat_curved <- lat_linear + offset * perp_y * len

  data.frame(lon = lon_curved, lat = lat_curved)
}

# Create curved paths for all edges
edges_with_paths <- edges_for_plot %>%
  rowwise() %>%
  mutate(
    path = list(generate_curved_path(lon_x, lat_x, lon_y, lat_y, n = 50, curvature = 0.15))
  ) %>%
  ungroup()

# Unnest the paths to get all points along each arc
edges_paths_long <- edges_with_paths %>%
  select(edge_pair, n_occurrences, path) %>%
  unnest(path) %>%
  group_by(edge_pair) %>%
  mutate(path_order = row_number()) %>%
  ungroup()

# Prepare node data for plotting
nodes_for_plot <- state_centroids %>%
  filter(state %in% states_involved) %>%
  left_join(REGION_DATA %>% select(state, bea_reg), by = "state")

# Add BEA regions to cam_map and filter out HI, AK, Mexico, and northern territories if not in network
cam_map_with_region <- cam_map %>%
  left_join(REGION_DATA %>% select(state, bea_reg), by = c("NAME_En" = "state")) %>%
  filter(
    !(NAME_En %in% c("Hawaii", "Alaska", "Mexico", "Yukon", "Northwest Territories", "Nunavut")) |
    NAME_En %in% states_involved
  )

# Create the map plot
p_map <- ggplot() +
  # Base map colored by region with black outlines
  geom_sf(data = cam_map_with_region, aes(fill = bea_reg), color = "black", linewidth = 0.5) +
  region_fill_scale() +
  # Add network edges as great circle arcs
  geom_path(data = edges_paths_long,
            aes(x = lon, y = lat, group = edge_pair,
                linewidth = n_occurrences, alpha = n_occurrences),
            color = "black") +
  # Add network nodes (circles with black outline)
  geom_point(data = nodes_for_plot,
             aes(x = lon, y = lat, fill = bea_reg),
             size = 5, shape = 21, color = "black", stroke = 1.5) +
  scale_linewidth_continuous(range = c(0.5, 3),
                             name = "# Time Points",
                             breaks = c(2, 5, 10, 15)) +
  scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
  guides(fill = "none") +  # Remove BEA region fill legend
  coord_sf(expand = FALSE) +  # Remove whitespace
  labs(title = "Conserved Significant Inter-Region Connections",
       subtitle = paste0("Connections appearing in ", min_occurrences, "+ time points (out of 20 quarters)\n",
                        "Network overlaid on BEA regions")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# Save map plot
fn_out_map <- paste0(fn_out_path, "network_conserved_", min_occurrences, "plus_map.jpg")
ggsave(fn_out_map, plot = p_map, width = 16, height = 10, dpi = 300)

cat("Map with network overlay saved to:", fn_out_map, "\n")

cat("\n=== Analysis Complete ===\n")
