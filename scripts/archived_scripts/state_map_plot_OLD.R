#Note this map code is custom built for CAM, given the differences of each region, code should probably be manually built for each

library(sf)
library(tidyverse)
library(patchwork)
source("scripts/color_schemes.R")

#Set as bounds for RR and transform to log for display purposes
b_list <- list(
  'Ontario' = c(0.1,0.2,1,5,15),
  'New York' = c(0.4,0.7,1,1.5,2),
  'British Columbia' = c(0.1,0.5,1,2,10),
  'Washington' = c(0.2,0.5,1,2,5),
  'Illinois' = c(0.4,0.7,1,1.5,3),
  'Georgia' = c(0.5,0.8,1,1.3,2),
  'Mexico' = c(0.05,0.1,1,10,20),
  'Texas' = c(0.4,0.7,1,1.5,2)
)
state_names <- names(b_list)

# Rotate an sf object around its centroid
rotate_geometry <- function(sf_obj, angle_degrees) {
  angle_radians <- angle_degrees * pi / 180
  geom <- sf::st_geometry(sf_obj)
  centroid <- sf::st_centroid(sf::st_union(geom))
  rotate_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians),
                            sin(angle_radians),  cos(angle_radians)),
                          nrow = 2)
  rotated_geom <- (geom - centroid) * rotate_matrix + centroid
  sf::st_set_geometry(sf_obj, rotated_geom)
}


df_RR_state <- read_tsv("results/CAM_1000/df_RR_by_state.tsv")
cam_map <- st_read("data/shp-files/cam-shp.gpkg")
st_bbox(cam_map)  #Should come back as meters (-640,000 for xmin)
st_crs(cam_map) <- 3857 #cam-shp is actually Web Mercator and encoded wrong by qGIS
cam_map <- cam_map %>%  
  st_transform(4326) #Transform back to CRS (WGS84)
st_bbox(cam_map) #Double check that it is actually DMS instead of meters

map_list <- list()
for(ORI in state_names){
  b <- b_list[ORI][[1]]
  
  #Subset to state to plot
  df_subset <- df_RR_state %>% 
    filter(x == ORI) %>%
    rename(state = y)
  df_subset$RR %>% quantile(c(0.2,0.95))
  
  
  df_map <- cam_map %>% 
    left_join(df_subset,by = join_by(NAME_En == state)) %>%
    mutate(fill_RR = fill_bound(RR,min(b),max(b)))
  
  df_main <- df_map %>% filter(NAME_En != "Hawaii") #Mainland
  
  
  # Extract Hawaii and rotate it 90 deg, then scale and translate geometry
  df_hi <- df_map %>%
    filter(NAME_En == "Hawaii") %>%
    rotate_geometry(45) %>%
    st_crop(xmin = -56,
            xmax = -47,
            ymin = -10,
            ymax = -5)
  
  ori_centroid <- df_main %>%
    filter(NAME_En == ORI) %>%
    st_point_on_surface() %>%
    st_transform(st_crs(df_main)) %>%
    st_coordinates() %>%
    as.data.frame() %>%
    mutate(label = "*")
  
  p_main <- ggplot(df_main) + 
    geom_sf(aes(fill = fill_RR), color = "black", size = 1) +
    RR_log_grad_manual(b) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(title = paste0(ORI)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  #Hawaii bbox with margins
  hi_bbox <- st_bbox(df_hi)
  hi_bbox["xmax"] <- hi_bbox["xmax"] + 1
  hi_bbox["xmin"] <- hi_bbox["xmin"] - (-0.5)
  hi_bbox["ymax"] <- hi_bbox["ymax"] + 0
  hi_bbox["ymin"] <- hi_bbox["ymin"] -1
  
  p_hi <- ggplot() + 
    geom_sf(data = st_as_sf(st_as_sfc(hi_bbox)), fill = "white", color = "black", size = 1.5) +
    geom_sf(data = df_hi, aes(fill = fill_RR), color = "black", size = 1) +
    RR_log_grad_manual(b) +
    theme_void() +
    theme(legend.position = "none")
  
  if(ORI == "Hawaii"){
    p_hi <- p_hi + 
      geom_text(data = ori_centroid, aes(X,Y,label=label),
                color = "black", size = 6, fontface = "bold")
  }else{
    p_main <- p_main + 
      geom_text(data = ori_centroid, aes(X,Y,label=label),
                color = "black", size = 6, fontface = "bold")
  }
  
  final_plot <- p_main +
    inset_element(p_hi, 
                  left = 0.05, bottom = 0.05, 
                  right = 0.25, top = 0.2)
  
  ggsave(filename = paste0("figs/CAM_1000/RR_maps/",ORI,".jpg"),
    plot = final_plot,
    height = 7,
    width = 7,
    units = "in",
    dpi = 192)
  map_list[[ORI]] <- final_plot
}

wrap_plots(map_list,ncol=2) +
  plot_annotation("Identical Sequence Enrichment") &
  theme(plot.title = element_text(
    hjust = 0.5,
    size = rel(1.5),
    face = "bold"
  ))
ggsave(filename="figs/CAM_1000/RR_maps/patched_maps.jpg",
       width = 15,
       height = 28,
       units = "in",
       dpi = 192)
