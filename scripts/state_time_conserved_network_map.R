library(tidyverse)
library(data.table)
library(argparse)
library(sf)

source("scripts/color_schemes.R")
source("scripts/cam_map.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--sd_threshold', type = 'double', default = 2.0, help = 'Number of standard deviations for significance threshold')
  parser$add_argument('--min_occurrences', type = 'integer', default = 4, help = 'Minimum number of time points a connection must appear in')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
sd_threshold <- args$sd_threshold
min_occurrences <- args$min_occurrences

select <- dplyr::select

# Read region data
REGION_DATA <- fread("data/us_states_regions.csv")

# Load CAM map using the existing function
cam_map <- prep_cam_map()

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

# Read the conserved connections data
fn_in <- paste0("results/", scenario, "/time_state/df_conserved_connections_", min_occurrences, "plus.tsv")
edge_counts <- fread(fn_in)

cat("\n=== Conserved Connection Map ===\n")
cat("Creating map overlay with", nrow(edge_counts), "conserved connections\n")

# Get unique states involved in conserved connections
states_involved <- unique(c(edge_counts$x, edge_counts$y))

# Prepare edge data for plotting
# Keep only unique edges (remove duplicates where x-y appears twice as x→y and y→x)
# Strategy: keep the edge where the northern state (higher latitude) is the source
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

# Create the plot
p <- ggplot() +
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

# Save plot
fn_out_path <- paste0("figs/", scenario, "/time_state_networks/")
fn_out <- paste0(fn_out_path, "network_conserved_", min_occurrences, "plus_map.jpg")
ggsave(fn_out, plot = p, width = 16, height = 10, dpi = 300)

cat("Map with network overlay saved to:", fn_out, "\n")
cat("\n=== Analysis Complete ===\n")
