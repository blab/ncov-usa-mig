# File: install_packages.R
# Description: Install required R packages using pak for the ncov-usa-mig analysis
# Date: 2025-10-02
# Updated: Complete package list from comprehensive scan of all scripts

# First, install pak if not already installed
if (!require("pak", quietly = TRUE)) {
  install.packages("pak", repos = "https://cloud.r-project.org")
}

# List of required packages
packages <- c(
  # Core data manipulation and tidyverse
  "tidyverse",      # Collection of data science packages (includes dplyr, ggplot2, readr, tidyr, purrr, stringr, tibble, forcats)
  "data.table",     # Fast data manipulation
  "dplyr",          # Data manipulation grammar
  "tidyr",          # Tidy messy data
  "readr",          # Fast file reading
  "purrr",          # Functional programming tools
  "stringr",        # String manipulation
  "lubridate",      # Date/time manipulation
  "magrittr",       # Pipe operators

  # Database connectivity
  "duckdb",         # DuckDB database connector
  "DBI",            # Database interface
  "dbplyr",         # Database operations with dplyr

  # Command-line interface
  "argparse",       # Command line argument parsing
  "cli",            # Terminal output formatting
  "rlang",          # Low-level R programming tools

  # Visualization - core
  "ggplot2",        # Grammar of graphics plotting
  "patchwork",      # Combine multiple plots
  "scales",         # Scale functions for visualization
  "RColorBrewer",   # Color palettes
  "viridis",        # Color-blind friendly palettes

  # Visualization - extensions
  "ggrepel",        # Non-overlapping text labels
  "ggpubr",         # Publication ready plots
  "ggsignif",       # Significance brackets for ggplot2
  "lemon",          # Axis and legend extensions for ggplot2
  "gridGraphics",   # Convert base graphics to grid

  # Network visualization
  "igraph",         # Network analysis and graph theory
  "ggraph",         # Grammar of graphics for network graphs

  # Spatial/geographic
  "sf",             # Simple features for spatial data
  "spdep",          # Spatial dependence analysis
  "tigris",         # Census geographic data
  "maps",           # Geographic map data
  "ggspatial",      # Spatial data visualization
  "usmap",          # US map plotting

  # Statistical modeling
  "broom",          # Tidying statistical model outputs
  "car",            # Companion to Applied Regression (VIF, bootstrapping)
  "lmtest",         # Linear model testing
  "sandwich",       # Robust covariance estimation
  "splines",        # Spline functions
  "mgcv",           # Generalized Additive Models
  "MASS",           # Robust linear models
  "quantreg",       # Quantile regression
  "lme4",           # Linear mixed-effects models
  "lmerTest",       # Testing linear mixed-effects models
  "rsample",        # For CV methods

  # Clustering and dimensionality reduction
  "cluster",        # Cluster analysis (agnes, diana, silhouette)
  "vegan",          # Community ecology (metaMDS, cmdscale)
  "factoextra",     # Visualization of clustering results
  "tsne",           # t-SNE dimensionality reduction
  "ggdendro",       # Dendrogram plotting with ggplot2
  "dendextend",     # Extended dendrogram functionality
  "mclust",         # Model-based clustering (Adjusted Rand Index)

  # Model evaluation
  "yardstick",      # Model performance metrics

  # Time series
  "zoo",            # Time series infrastructure

  # Utilities
  "jsonlite",       # JSON parsing and generation
  "magick",         # Image processing
  "tictoc",         # Timing code execution
  "doParallel",     # Parallel processing backend

  # 3D visualization
  "plotly",         # Interactive 3D plots
  "htmlwidgets",    # Save interactive HTML widgets
  "scatterplot3d"   # Static 3D scatter plots
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