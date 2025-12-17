# --- deps ---
library(sf)
library(tidyverse)
library(patchwork)     # for inset_element()
# source("scripts/color_schemes.R") # keep this if you use your custom scales

# --- duct-tape geometry helper: rotate around centroid (kept as-is) ---
rotate_geometry <- function(sf_obj, angle_degrees) {
  angle_radians <- angle_degrees * pi / 180
  geom <- sf::st_geometry(sf_obj)
  centroid <- sf::st_centroid(sf::st_union(geom))
  R <- matrix(c(cos(angle_radians), -sin(angle_radians),
                sin(angle_radians),  cos(angle_radians)), nrow = 2)
  rotated_geom <- (geom - centroid) * R + centroid
  sf::st_set_geometry(sf_obj, rotated_geom)
}

# --- load CAM map and fix mislabeled CRS, then reproject to WGS84 ---
prep_cam_map <- function(gpkg_path = "data/shp-files/cam-shp.gpkg") {
  cam_map <- st_read(gpkg_path, quiet = TRUE)
  # "duct tape": file is actually Web Mercator (EPSG:3857) but encoded wrong
  st_crs(cam_map) <- 3857          # relabel only (no reprojection)
  cam_map <- st_transform(cam_map, 4326)  # proper reprojection to WGS84
  cam_map
}

# --- main reusable choropleth with HI inset ---
# data: tibble with state names and a value column to map
# state_col: column in `data` that matches cam_map$NAME_En
# value_col: column in `data` to color by
# fill_mapper: function(raw_vector) -> mapped fill (default identity)
# scale_fun: a function with no args that returns a ggplot2 scale (e.g., RR_log_grad_manual(b))
# ori: optional state name to mark with an asterisk
plot_cam_choropleth <- function(
    cam_map,
    data,
    state_col,
    value_col,
    fill_mapper = identity, 
    scale_fun = NULL,
    title = NULL,
    ori = NULL,
    hi_angle = 45,
    inset_loc = list(left = 0.05, bottom = 0.05, right = 0.25, top = 0.2),
    hi_crop = c(xmin = -56, xmax = -47, ymin = -10, ymax = -5),
    bbox_margin = c(xmin = -0.5, xmax = 1, ymin = -1, ymax = 0),
    line_size = 1, #For state borderlines
    box_size = 1.5, #For box/frame around HI
    bottom_legend = FALSE,
    theme_override = NULL
) {
  state_sym <- rlang::ensym(state_col)
  value_sym <- rlang::ensym(value_col)
  
  # join map with data and create mapped fill column
  df_map <- cam_map %>%
    left_join(data, by = join_by(NAME_En == !!state_sym)) %>%
    mutate(.fill_value = fill_mapper(!!value_sym))
  
  # split main vs HI
  df_main <- df_map %>% filter(NAME_En != "Hawaii")
  df_hi <- df_map %>%
    filter(NAME_En == "Hawaii") %>%
    rotate_geometry(45) %>%
    st_crop(xmin = -56,
            xmax = -47,
            ymin = -10,
            ymax = -5)
  
  # optional asterisk for origin
  ori_centroid <- NULL
  if (!is.null(ori)) {
    ori_centroid <- df_map %>%
      filter(NAME_En == ori) %>%
      st_point_on_surface() %>%
      st_coordinates() %>%
      as.data.frame() %>%
      setNames(c("X", "Y")) %>%
      mutate(label = "*")
  }
  
  # main plot
  p_main <- ggplot(df_main) +
    geom_sf(aes(fill = .fill_value), color = "black", size = line_size) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  # Apply theme override if provided
  if (!is.null(theme_override)) {
    p_main <- p_main + theme_override
  }
  
  # add custom fill scale if provided
  if (is.function(scale_fun)) {
    p_main <- p_main + scale_fun()
  }
  
  if (!is.null(title)) {
    p_main <- p_main + labs(title = title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  if(bottom_legend) {
    p_main <- p_main +
      theme(legend.position = "bottom") +
      theme(legend.text = element_text(angle = 45, hjust = 1))
  } else{
    p_main <- p_main +
      theme(legend.position = "right")
  }
  
  # place asterisk either on main or on HI depending on ori
  add_star_to_main <- !is.null(ori) && ori != "Hawaii"
  
  if (add_star_to_main && !is.null(ori_centroid)) {
    p_main <- p_main +
      geom_text(
        data = ori_centroid,
        aes(X, Y, label = label),
        inherit.aes = FALSE,
        color = "black", size = 6, fontface = "bold"
      )
  }
  
  # HI inset frame with margins
  hi_bbox <- st_bbox(df_hi)
  hi_bbox["xmax"] <- hi_bbox["xmax"] + bbox_margin["xmax"]
  hi_bbox["xmin"] <- hi_bbox["xmin"] - (-bbox_margin["xmin"])
  hi_bbox["ymax"] <- hi_bbox["ymax"] + bbox_margin["ymax"]
  hi_bbox["ymin"] <- hi_bbox["ymin"] - (-bbox_margin["ymin"])
  hi_frame <- st_as_sf(st_as_sfc(hi_bbox))
  
  p_hi <- ggplot() +
    geom_sf(data = hi_frame, fill = "white", color = "black", size = box_size) +
    geom_sf(data = df_hi, aes(fill = .fill_value), color = "black", size = line_size) +
    theme_void() +
    theme(legend.position = "none") 
  
  if (is.function(scale_fun)) {
    p_hi <- p_hi + scale_fun()
  }
  
  if (!add_star_to_main && !is.null(ori_centroid)) {
    p_hi <- p_hi +
      geom_text(
        data = ori_centroid,
        aes(X, Y, label = label),
        inherit.aes = FALSE,
        color = "black", size = 6, fontface = "bold"
      )
  }
  
  # compose with inset
  p_main +
    inset_element(
      p_hi,
      left = inset_loc$left, bottom = inset_loc$bottom,
      right = inset_loc$right, top = inset_loc$top,
      ignore_tag = TRUE
    )
}