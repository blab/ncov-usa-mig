#File: make_state_nb_dist.R
#Author(s): Amin Bemanian
#Date: 8/19/24
#Description: R script that obtains state/province shapefiles from cam-shp.gpkg and makes a
#             neighbor (adjacency) and euclidean distance matrix using spdep/sf.
#Does not need to be re-run with scenarios after initial run.

library(spdep)
library(sf)
library(dplyr)
library(data.table)

# Turn off S2 spherical geometry to avoid geometry validation issues
sf_use_s2(FALSE)

MAX_LAG <- 12 #Maximum higher order neighbors, California and Maine have the farthest queen adjacency distance at 11, may need to change if you use rook adjacency

# Read CAM shapefile (USA + Canada + Mexico)
message(Sys.time(), " - Reading CAM shapefile...")
cam_shp <- st_read("data/shp-files/cam-shp.gpkg", layer = "union", quiet = TRUE)
message(Sys.time(), " - Finished reading shapefile")

# Check number of polygons
n_polygons <- nrow(cam_shp)
message(Sys.time(), " - Number of polygons: ", n_polygons)
message("Expected: 65 polygons (states/provinces)")

# Check geometry validity
message(Sys.time(), " - Checking geometry validity...")
valid_geom <- st_is_valid(cam_shp)
n_invalid <- sum(!valid_geom)
message("Invalid geometries: ", n_invalid)
if(n_invalid > 0){
  message("Invalid polygon names:")
  print(cam_shp$NAME_En[!valid_geom])
}

# Simplify geometries to speed up poly2nb if needed
# NOTE: preserveTopology = FALSE is used for performance. With complex Arctic geometries
# (Nunavut, Alaska), topology-preserving simplification hangs indefinitely.
# The poly2nb() snap parameter handles any small gaps created by non-topological simplification.
# Tolerance is in DECIMAL DEGREES for lon/lat data: 0.1 degrees ≈ 11 km at equator (~8 km at 49°N)
message(Sys.time(), " - Simplifying geometries (tolerance = 0.1 degrees ≈ 11km)...")
cam_shp <- st_simplify(cam_shp, preserveTopology = FALSE, dTolerance = 1)
message(Sys.time(), " - Geometry simplification complete")

# Create neighbor list using queen adjacency
message(Sys.time(), " - Creating neighbor list using queen adjacency...")
STATE_NAMES <- cam_shp$NAME_En
nb_states <- poly2nb(cam_shp, row.names = STATE_NAMES)

# Manual correction: Idaho-BC border is lost during simplification
# The Idaho-BC border is very small (~45 miles) and gets removed even with conservative tolerances
# They are verified to touch in the original unsimplified shapefile
message(Sys.time(), " - Applying manual corrections for borders lost during simplification...")
idaho_idx <- which(STATE_NAMES == "Idaho")
bc_idx <- which(STATE_NAMES == "British Columbia")
if(length(idaho_idx) == 1 && length(bc_idx) == 1) {
  # Add BC as a neighbor of Idaho if not already present
  if(!(bc_idx %in% nb_states[[idaho_idx]])) {
    nb_states[[idaho_idx]] <- c(nb_states[[idaho_idx]], bc_idx)
    message("  - Added British Columbia as neighbor of Idaho")
  }
  # Add Idaho as a neighbor of BC if not already present
  if(!(idaho_idx %in% nb_states[[bc_idx]])) {
    nb_states[[bc_idx]] <- c(nb_states[[bc_idx]], idaho_idx)
    message("  - Added Idaho as neighbor of British Columbia")
  }
}
message(Sys.time(), " - Manual corrections complete")

# Compute lagged neighbors
message(Sys.time(), " - Computing lagged neighbors (max lag = ", MAX_LAG, ")...")
nb_lagged_states <- nblag(nb_states, MAX_LAG)
message(Sys.time(), " - Finished neighbor list computation (", length(STATE_NAMES), " states/provinces)")

# Switch from sparse lists to a long matrix
message(Sys.time(), " - Creating state pair matrix...")
nb_matrix <- expand.grid(STATE_NAMES, STATE_NAMES, stringsAsFactors = FALSE) %>% as.data.table
colnames(nb_matrix) <- c("state_x", "state_y")

nb_matrix$nb_dist <- as.numeric(-1) #Temporarily set default distance to -1 to initialize

# Calculate neighbor distances using vectorized approach
# Pre-allocate matrix with correct dimensions
message(Sys.time(), " - Calculating neighbor distances...")
n_states <- length(STATE_NAMES)
nb_dist_matrix <- matrix(-1, nrow = n_states, ncol = n_states)
diag(nb_dist_matrix) <- 0  # Set diagonal to 0 (self-distance)

# Populate neighbor distances
for(l in 1:MAX_LAG){
  nb_level <- nb_lagged_states[[l]]
  for(i in 1:n_states){
    nb_i <- nb_level[[i]]
    if(length(nb_i) > 0){
      nb_dist_matrix[i, nb_i] <- l
    }
  }
}
message(Sys.time(), " - Finished neighbor distance calculation")

# Calculate euclidean distances between state/province centroids
message(Sys.time(), " - Calculating state centroids...")
state_centroids <- st_centroid(cam_shp)
message(Sys.time(), " - Calculating euclidean distances between centroids...")
state_dists <- st_distance(state_centroids) * 1E-3 #Convert to km
message(Sys.time(), " - Finished euclidean distance calculation")

# Convert matrices to long format using vectorized operations
message(Sys.time(), " - Converting matrices to long format...")
nb_dist_vec <- as.vector(nb_dist_matrix)
nb_dist_vec[nb_dist_vec == -1] <- NA  # Convert -1 to NA
euclid_dist_vec <- as.vector(state_dists)

# Assign to nb_matrix using vectorized operations
nb_matrix$nb_dist <- nb_dist_vec
nb_matrix$euclid_dist <- euclid_dist_vec

# Order and save
message(Sys.time(), " - Ordering and writing output file...")
nb_matrix <- nb_matrix[order(nb_matrix$state_x, nb_matrix$state_y),]
readr::write_tsv(nb_matrix, "data/nb_dist_states.tsv")

message("Neighbor distance matrix saved to data/nb_dist_states.tsv")
message("Total state/province pairs: ", nrow(nb_matrix))
