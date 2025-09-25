source("scripts/cam_map.R")
source("scripts/color_schemes.R")

# Load shapefile and RR file
cam_map <- prep_cam_map("data/shp-files/cam-shp.gpkg")
df_RR_state <- read_tsv("results/CAM_1000/df_RR_by_state.tsv")

# List of states/provinces to use as origin points
b_list <- list(
  'Ontario' = c(0.1,0.2,1,5,15),
  'New York' = c(0.4,0.7,1,1.5,2),
  'British Columbia' = c(0.1,0.5,1,2,10),
  'Washington' = c(0.2,0.5,1,2,5),
  'Illinois' = c(0.4,0.7,1,1.5,3),
  'Georgia' = c(0.1,0.5,1,2,10),
  'Mexico' = c(0.02,0.1,1,10,50),
  'Texas' = c(0.4,0.7,1,1.5,2)
)

map_list <- vector("list", length(b_list))
names(map_list) <- names(b_list)

# Loop over each and make a map
for (ORI in names(b_list)) {
  b <- b_list[[ORI]]
  
  df_subset <- df_RR_state %>%
    filter(x == ORI) %>%
    rename(state = y)
  
  fill_mapper <- function(rr) fill_bound(rr, min(b), max(b))
  scale_fun <- function() RR_log_grad_manual(b)
  
  p <- plot_cam_choropleth(
    cam_map = cam_map,
    data = df_subset,
    state_col = state,
    value_col = RR,
    fill_mapper = fill_mapper,
    scale_fun = scale_fun,
    title = ORI,
    ori = ORI
  )
  
  ggsave(paste0("figs/CAM_1000/RR_maps/", ORI, ".jpg"),
         plot = p, width = 7, height = 7, units = "in", dpi = 192)
  
  map_list[[ORI]] <- p
}

# Combine into a single polot
wrap_plots(map_list, ncol = 2) +
  plot_annotation(
  title="Identical Sequence Enrichment",
  theme = theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.5), face = "bold")
  ))
ggsave(filename="figs/CAM_1000/RR_maps/patched_maps.jpg",
       width = 15,
       height = 28,
       units = "in",
       dpi = 192)
 