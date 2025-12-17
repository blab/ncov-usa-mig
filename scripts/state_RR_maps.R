source("scripts/cam_map.R")
source("scripts/color_schemes.R")

# Load shapefile and RR file
cam_map <- prep_cam_map("data/shp-files/cam-shp.gpkg")
df_RR_state <- read_tsv("results/CAM_1000/df_RR_by_state.tsv")

# List of states/provinces to use as origin points
b_list <- list( #Bounds
  'Ontario' = c(0.1,0.2,1,5,15),
  'New York' = c(0.4,0.7,1,1.5,2),
  #'California' = c(0.2,0.5,1,2,5),
  'Mexico' = c(0.02,0.1,1,10,50)
)

place_abbrv <- list( #Abbreviatons
  'Ontario' = 'ON',
  'New York' = 'NY',
  #'California' = 'CA',
  'Mexico' = 'MX'
)

map_list <- vector("list", length(b_list))
names(map_list) <- names(b_list)

# Loop over each and make a map
for (ORI in names(b_list)) {
  b <- b_list[[ORI]]
  abrv <- place_abbrv[[ORI]]
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
    title = abrv,
    ori = ORI,
    line_size = 0.1,
    box_size = 0.2
  )
  
  ggsave(paste0("figs/CAM_1000/RR_maps/", ORI, ".jpg"),
         plot = p, width = 3, height = 3, units = "in", dpi = 300)
  
  map_list[[ORI]] <- p
}

# Combine into a single polot
wrap_plots(map_list, nrow = 1) +
  plot_annotation(
  title="Identical Sequence Relative Risk",
  theme = theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.5), face = "bold")
  ))
ggsave(filename="figs/CAM_1000/RR_maps/patched_maps.jpg",
       width = 15,
       height = 4,
       units = "in",
       dpi = 300)
 