# File: plot_region_map.R
# Author(s): ChatGPT, based on color_schemes.R by Amin Bemanian
# Date: 2025-09-26
# Description: Create a map visualization of US BEA regions using defined color schemes

library(sf)
library(tidyverse)
library(patchwork)
library(argparse)
source("scripts/color_schemes.R")
source("scripts/cam_map.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# Read BEA region assignments
regions_data <- read_csv("data/us_states_regions.csv", show_col_types = FALSE) %>%
  select(state, region = bea_reg) %>%
  # Ensure state names match cam_map
  mutate(
    state = case_when(
      state == "District of Columbia" ~ "Washington DC",
      TRUE ~ state
    ),
    # Order the legend to match REGION_SCALE's color order instead of alphabetical
    region = factor(region, levels = names(REGION_SCALE))
  )

# Prepare the map
cam_map <- prep_cam_map()

# Create the plot using cam_map.R functionality
p <- plot_cam_choropleth(
  cam_map = cam_map,
  data = regions_data,
  state_col = state,
  value_col = region,
  fill_mapper = identity,
  scale_fun = region_fill_scale,  # Use the predefined scale
  title = NULL,
  bottom_legend = TRUE,
  legend_angle = 0,  # region names read better horizontal than the default 45deg
  guide_override = guides(fill = guide_legend(nrow = 3, byrow = TRUE)),
  theme_override = theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )
)

# Create output directory if it doesn't exist
dir.create(paste0("figs/", scenario), recursive = TRUE, showWarnings = FALSE)

# Save the plot (wide canvas so the 11-item legend has room to lay out
# without clipping; the stitch script rescales based on aspect ratio, so
# this native size doesn't constrain the final printed size)
ggsave(
  paste0("figs/", scenario, "/bea_region_map.png"),
  p,
  width = 12.5,
  height = 8.5,
  units = "in",
  dpi = 300
)