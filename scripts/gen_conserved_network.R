library(tidyverse)
library(data.table)
library(argparse)

source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  parser$add_argument('--percentile_threshold', type = 'double', default = 0.98, help = 'Percentile threshold for significance (default is 0.98)')
  parser$add_argument('--min_occurrences', type = 'integer', default = 5, help = 'Minimum number of time points a connection must appear in')
  parser$add_argument('--gen_viz', action = 'store_true', default = FALSE, help = 'Generate network visualizations (off by default)')
  return(parser$parse_args())
}

args <- collect_args()
scenario          <- args$scenario
percentile_threshold <- args$percentile_threshold
min_occurrences   <- args$min_occurrences
gen_viz           <- args$gen_viz
label_threshold   <- round(100 * percentile_threshold)

select <- dplyr::select

REGION_DATA <- fread("data/us_states_regions.csv")

# Read the significant connections
fn_in <- paste0("results/", scenario, "/time_state/df_significant_connections_", label_threshold, "_percentile.tsv")
sig_connections <- fread(fn_in) %>%
  mutate(edge_pair = paste(pmin(x, y), pmax(x, y), sep = "_"))

# Count how many times each edge appears across all time points
edge_counts <- sig_connections %>%
  group_by(edge_pair, x, y, region_x, region_y, country_x, country_y) %>%
  summarise(
    n_occurrences = n(),
    first_date    = min(date),
    last_date     = max(date),
    mean_nRR      = mean(nRR, na.rm = TRUE),
    median_nRR    = median(nRR, na.rm = TRUE),
    max_nRR       = max(nRR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_occurrences >= min_occurrences) %>%
  filter(region_x != region_y) %>%
  filter(x < y) %>%
  arrange(desc(n_occurrences), desc(mean_nRR))

cat("\n=== Conserved Connection Analysis ===\n")
cat("Percentile threshold:", label_threshold, "%\n")
cat("Connections appearing in", min_occurrences, "or more time points:", nrow(edge_counts), "\n")
cat("\nTop 10 most conserved connections:\n")
print(head(edge_counts %>% select(x, y, n_occurrences, first_date, last_date, mean_nRR), 10))

# Save the conserved edges data
fn_out_path <- paste0("results/", scenario, "/time_state/")
fn_out_data <- paste0(fn_out_path, "df_conserved_connections_", label_threshold, "pct_", min_occurrences, "plus.tsv")
write_tsv(edge_counts, fn_out_data)
cat("Conserved connections data saved to:", fn_out_data, "\n")

if (gen_viz) {
  library(igraph)
  library(ggraph)
  library(sf)
  library(patchwork)

  # Load CAM shapefile and extract centroids
  cam_map <- st_read("data/shp-files/cam-shp.gpkg", quiet = TRUE)
  st_crs(cam_map) <- 3857
  cam_map <- st_transform(cam_map, 4326)

  state_centroids <- cam_map %>%
    st_centroid() %>%
    mutate(
      lon = st_coordinates(.)[, 1],
      lat = st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry() %>%
    as_tibble() %>%
    select(state = NAME_En, lon, lat)

  # ============================================================================
  # IGRAPH NETWORK VISUALIZATION
  # ============================================================================

  conserved_edges <- edge_counts %>%
    select(x, y, n_occurrences, mean_nRR, max_nRR)

  g <- graph_from_data_frame(conserved_edges, directed = FALSE)

  node_names <- V(g)$name
  node_data  <- REGION_DATA %>%
    filter(state %in% node_names) %>%
    select(state, bea_reg, country)
  node_data <- node_data[match(node_names, node_data$state), ]

  V(g)$region  <- node_data$bea_reg
  V(g)$country <- node_data$country

  node_coords <- state_centroids %>%
    filter(state %in% node_names)
  node_coords <- node_coords[match(node_names, node_coords$state), ]

  set.seed(42)
  jitter_amount <- 0.3
  node_coords <- node_coords %>%
    mutate(
      lon = lon + runif(n(), -jitter_amount, jitter_amount),
      lat = lat + runif(n(), -jitter_amount, jitter_amount)
    )

  layout_matrix <- as.matrix(node_coords[, c("lon", "lat")])

  E(g)$weight   <- conserved_edges$n_occurrences
  E(g)$max_nRR  <- conserved_edges$max_nRR

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
         subtitle = paste0("Connections above ", label_threshold, "th percentile in ", min_occurrences, "+ time points\n",
                          "Total conserved connections: ", nrow(conserved_edges), " | States involved: ", vcount(g))) +
    theme_graph() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.position = "right"
    )

  fn_viz_path <- paste0("figs/", scenario, "/time_state_networks/")
  fn_out_net  <- paste0(fn_viz_path, "network_conserved_", label_threshold, "pct_", min_occurrences, "plus.jpg")
  ggsave(fn_out_net, plot = p, width = 14, height = 9, dpi = 300, create.dir = TRUE)
  cat("\nConserved network plot saved to:", fn_out_net, "\n")

  # ============================================================================
  # GEOGRAPHIC MAP VISUALIZATION
  # ============================================================================
  cat("\n=== Creating Geographic Map Overlay ===\n")

  # Drop Canadian provinces from the geographic map panel (left panel)
  canada_states <- REGION_DATA %>% filter(country == "Canada") %>% pull(state)
  edge_counts_map <- edge_counts %>%
    filter(!(x %in% canada_states), !(y %in% canada_states))

  states_involved <- unique(c(edge_counts_map$x, edge_counts_map$y))

  generate_curved_path <- function(lon1, lat1, lon2, lat2, n = 20, curvature = 0.15) {
    t <- seq(0, 1, length.out = n)
    lon_linear <- lon1 + t * (lon2 - lon1)
    lat_linear <- lat1 + t * (lat2 - lat1)
    offset <- curvature * sin(t * pi)
    dx <- lon2 - lon1
    dy <- lat2 - lat1
    len <- sqrt(dx^2 + dy^2)
    perp_x <- -dy / len
    perp_y <-  dx / len
    lon_curved <- lon_linear + offset * perp_x * len
    lat_curved <- lat_linear + offset * perp_y * len
    data.frame(lon = lon_curved, lat = lat_curved)
  }

  edges_for_plot <- edge_counts_map %>%
    left_join(state_centroids, by = c("x" = "state")) %>%
    rename(lon_x = lon, lat_x = lat) %>%
    left_join(state_centroids, by = c("y" = "state")) %>%
    rename(lon_y = lon, lat_y = lat) %>%
    group_by(edge_pair) %>%
    slice_min(lon_x, n = 1, with_ties = FALSE) %>%
    ungroup()

  edges_paths_long <- edges_for_plot %>%
    rowwise() %>%
    mutate(path = list(generate_curved_path(lon_x, lat_x, lon_y, lat_y, n = 20, curvature = 0.15))) %>%
    ungroup() %>%
    select(edge_pair, n_occurrences, path) %>%
    unnest(path) %>%
    group_by(edge_pair) %>%
    mutate(path_order = row_number()) %>%
    ungroup()

  nodes_for_plot <- state_centroids %>%
    filter(state %in% states_involved) %>%
    left_join(REGION_DATA %>% select(state, bea_reg), by = "state")

  cam_map_with_region <- cam_map %>%
    left_join(REGION_DATA %>% select(state, bea_reg), by = c("NAME_En" = "state")) %>%
    filter(
      !(NAME_En %in% canada_states),
      !(NAME_En %in% c("Hawaii", "Alaska", "Mexico", "Yukon", "Northwest Territories", "Nunavut")) |
      NAME_En %in% states_involved
    )

  p_map <- ggplot() +
    geom_sf(data = cam_map_with_region, aes(fill = bea_reg), color = "black", linewidth = 0.5) +
    region_fill_scale() +
    geom_path(data = edges_paths_long,
              aes(x = lon, y = lat, group = edge_pair,
                  linewidth = n_occurrences, alpha = n_occurrences),
              color = "black") +
    geom_point(data = nodes_for_plot,
               aes(x = lon, y = lat, fill = bea_reg),
               size = 2, shape = 21, color = "black", stroke = 0.5) +
    scale_linewidth_continuous(range = c(0.5, 3),
                               name = "# Time Points",
                               breaks = c(2, 5, 10, 15)) +
    scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
    guides(fill = "none") +
    coord_sf(expand = FALSE) +
    labs(title = paste0("All Quarters (n = ", nrow(conserved_edges), ")")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )

  fn_out_map <- paste0(fn_viz_path, "network_conserved_", label_threshold, "pct_", min_occurrences, "plus_map.jpg")
  ggsave(fn_out_map, plot = p_map, width = 16, height = 10, dpi = 300)
  cat("Map with network overlay saved to:", fn_out_map, "\n")

  # ============================================================================
  # QUARTERLY NETWORK VISUALIZATION (combined 2x2 facet)
  # ============================================================================
  cat("\n=== Creating Combined Quarterly Network Map ===\n")

  EXCLUDE_STATES <- c("Hawaii", "Alaska", "Mexico", "Canada",
                      "Alberta", "British Columbia", "Manitoba", "New Brunswick",
                      "Newfoundland and Labrador", "Northwest Territories", "Nova Scotia",
                      "Nunavut", "Ontario", "Prince Edward Island", "Quebec",
                      "Saskatchewan", "Yukon")

  sig_connections_q <- sig_connections %>%
    filter(!(x %in% EXCLUDE_STATES | y %in% EXCLUDE_STATES)) %>%
    mutate(quarter = paste0("Q", lubridate::quarter(date)))

  quarter_month_labels <- c("Q1" = "Jan-Mar", "Q2" = "Apr-Jun",
                            "Q3" = "Jul-Sep", "Q4" = "Oct-Dec")
  min_occ_quarterly <- 2

  cam_map_conus <- cam_map %>%
    left_join(REGION_DATA %>% select(state, bea_reg), by = c("NAME_En" = "state")) %>%
    filter(!(NAME_En %in% EXCLUDE_STATES))

  quarterly_edges <- list()
  quarterly_nodes <- list()

  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    edge_counts_q <- sig_connections_q %>%
      filter(quarter == q) %>%
      group_by(edge_pair, x, y, region_x, region_y) %>%
      summarise(n_occurrences = n(), .groups = "drop") %>%
      filter(n_occurrences >= min_occ_quarterly) %>%
      filter(region_x != region_y)

    cat("  ", q, ":", nrow(edge_counts_q), "connections\n")
    if (nrow(edge_counts_q) == 0) next

    edges_for_plot_q <- edge_counts_q %>%
      left_join(state_centroids, by = c("x" = "state")) %>%
      rename(lon_x = lon, lat_x = lat) %>%
      left_join(state_centroids, by = c("y" = "state")) %>%
      rename(lon_y = lon, lat_y = lat) %>%
      group_by(edge_pair) %>%
      slice_min(lon_x, n = 1, with_ties = FALSE) %>%
      ungroup()

    edges_paths_long_q <- edges_for_plot_q %>%
      rowwise() %>%
      mutate(path = list(generate_curved_path(lon_x, lat_x, lon_y, lat_y,
                                              n = 20, curvature = 0.15))) %>%
      ungroup() %>%
      select(edge_pair, n_occurrences, path) %>%
      unnest(path) %>%
      mutate(quarter = q)

    states_involved_q <- unique(c(edge_counts_q$x, edge_counts_q$y))
    nodes_q <- state_centroids %>%
      filter(state %in% states_involved_q) %>%
      left_join(REGION_DATA %>% select(state, bea_reg), by = "state") %>%
      mutate(quarter = q)

    quarterly_edges[[q]] <- edges_paths_long_q
    quarterly_nodes[[q]] <- nodes_q
  }

  all_q_edges <- bind_rows(quarterly_edges) %>%
    mutate(n_occurrences = pmin(n_occurrences, 4))

  quarter_counts <- all_q_edges %>%
    distinct(quarter, edge_pair) %>%
    count(quarter)
  quarter_facet_labels <- setNames(
    paste0(quarter_month_labels[quarter_counts$quarter],
           " (n = ", quarter_counts$n, ")"),
    quarter_counts$quarter
  )
  for (q in c("Q1", "Q2", "Q3", "Q4")) {
    if (!q %in% names(quarter_facet_labels))
      quarter_facet_labels[[q]] <- paste0(quarter_month_labels[q], " (n = 0)")
  }
  quarter_facet_labels <- quarter_facet_labels[c("Q1", "Q2", "Q3", "Q4")]

  all_q_edges <- all_q_edges %>%
    mutate(quarter = factor(quarter, levels = c("Q1","Q2","Q3","Q4"),
                            labels = quarter_facet_labels))
  all_q_nodes <- bind_rows(quarterly_nodes) %>%
    mutate(quarter = factor(quarter, levels = c("Q1","Q2","Q3","Q4"),
                            labels = quarter_facet_labels))

  p_quarters <- ggplot() +
    geom_sf(data = cam_map_conus, aes(fill = bea_reg),
            color = "black", linewidth = 0.2) +
    region_fill_scale() +
    geom_path(data = all_q_edges,
              aes(x = lon, y = lat,
                  group = interaction(quarter, edge_pair),
                  linewidth = n_occurrences,
                  alpha = n_occurrences),
              color = "black") +
    geom_point(data = all_q_nodes,
               aes(x = lon, y = lat, fill = bea_reg),
               size = 1.2, shape = 21, color = "black", stroke = 0.3) +
    scale_linewidth_continuous(range = c(0.3, 3),
                               name = "# Years",
                               breaks = c(2, 3, 4),
                               limits = c(2, 4)) +
    scale_alpha_continuous(range = c(0.5, 0.9), guide = "none") +
    guides(fill = "none") +
    facet_wrap(~ quarter, nrow = 2, ncol = 2) +
    coord_sf(expand = FALSE) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(2, 2, 2, 2),
      panel.spacing = unit(2, "pt"),
      strip.text = element_text(size = 10, face = "bold")
    )

  fn_out_q <- paste0(fn_viz_path, "network_conserved_", label_threshold, "pct_quarterly_map.jpg")
  ggsave(fn_out_q, plot = p_quarters, width = 8, height = 8, dpi = 300)
  cat("Quarterly faceted map saved to:", fn_out_q, "\n")

  # ============================================================================
  # COMBINED FIGURE
  # ============================================================================
  combined <- p_map + p_quarters + plot_layout(ncol = 2) +
    plot_annotation(
      title = "Conserved Inter-Region Connections",
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
    )
  fn_out_combined     <- paste0(fn_viz_path, "network_conserved_combined.jpg")
  fn_out_combined_svg <- paste0(fn_viz_path, "network_conserved_combined.svg")
  ggsave(fn_out_combined,     plot = combined, width = 7, height = 4, dpi = 300)
  ggsave(fn_out_combined_svg, plot = combined, width = 7, height = 4)
  cat("Combined figure saved to:", fn_out_combined, "\n")
  cat("Combined SVG saved to:",    fn_out_combined_svg, "\n")
}

cat("\n=== Analysis Complete ===\n")
