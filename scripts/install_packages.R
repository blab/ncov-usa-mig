# File: install_packages.R
# Description: Install required R packages using pak for the ncov-usa-mig analysis
# Date: 2025-09-25

# First, install pak if not already installed
if (!require("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cloud.r-project.org")
}

# List of required packages
packages <- c(
  # Core data manipulation and visualization
  "tidyverse",      # Collection of data science packages
  "data.table",     # Fast data manipulation
  "dplyr",          # Data manipulation grammar
  "dbplyr",         # Database operations with dplyr
  "readr",          # Fast file reading
  
  # Database connectivity
  "duckdb",        # DuckDB database connector
  "DBI",           # Database interface
  
  # Visualization
  "ggplot2",       # Grammar of graphics plotting
  "scales",        # Scale functions for visualization
  "viridis",       # Color palettes
  "patchwork",     # Combine plots
  "ggpubr",        # Publication ready plots
  "magick",        # Image manipulation
  
  # Statistical analysis
  "splines",       # Regression spline functions
  "broom",         # Convert statistical objects to tidy format
  "sandwich",      # Robust covariance matrix estimators
  "lmtest",        # Testing Linear Regression Models
  "car",           # Companion to Applied Regression
  "mclust",        # Clustering algorithms
  "spdep",         # Spatial dependence
  "yardstick",     # Model metric computations
  
  # Spatial data handling
  "maps",          # Map visualization
  "sf",            # Simple features for spatial data
  "ape",           # Analysis of Phylogenetics and Evolution
  
  # Utility and performance
  "tictoc",        # Timing functions
  "argparse",      # Command line argument parsing
  "doParallel",    # Parallel backend
  
  # Additional tidyverse components
  "tidyr",         # Tidy messy data
  "purrr"          # Functional programming tools
)

# Install packages using pak
message("Installing required packages...")
pak::pkg_install(packages)

# Verify installation
installed <- rownames(installed.packages())
missing <- packages[!packages %in% installed]

if (length(missing) > 0) {
  warning("The following packages could not be installed: ", 
          paste(missing, collapse = ", "))
} else {
  message("All packages successfully installed!")
}

# Print package versions for reproducibility
message("\nInstalled package versions:")
for (pkg in packages) {
  if (pkg %in% installed) {
    version <- packageVersion(pkg)
    message(sprintf("- %s: %s", pkg, version))
  }
}