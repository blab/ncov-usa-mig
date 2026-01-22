library(sf)
library(tidyverse)
library(patchwork)
library(duckdb)

source("scripts/cam_map.R")

# Connect to database
con <- DBI::dbConnect(duckdb(), "db_files/db_CAM_1000.duckdb")

# Get reporting rates by division
df_reporting <- tbl(con, "metadata") %>%
  group_by(division) %>%
  summarise(
    total = n(),
    valid_age = sum(ifelse(!is.na(age_adj), 1, 0), na.rm = TRUE),
    valid_sex = sum(ifelse(sex != "", 1, 0), na.rm = TRUE)
  ) %>%
  collect() %>%
  mutate(
    pct_age = valid_age / total * 100,
    pct_sex = valid_sex / total * 100
  )

DBI::dbDisconnect(con, shutdown = TRUE)

# Load map
cam_map <- prep_cam_map()

# Create color scale for percentage
pct_scale <- function() {
  scale_fill_viridis_c(
    name = "% Valid",
    limits = c(0, 100),
    option = "plasma"
  )
}

# Age reporting map
p_age <- plot_cam_choropleth(
  cam_map = cam_map,
  data = df_reporting,
  state_col = division,
  value_col = pct_age,
  scale_fun = pct_scale,
  title = "Age Data Reporting Rate"
)

# Sex reporting map
p_sex <- plot_cam_choropleth(
  cam_map = cam_map,
  data = df_reporting,
  state_col = division,
  value_col = pct_sex,
  scale_fun = pct_scale,
  title = "Sex Data Reporting Rate"
)

# Combine maps
p_combined <- p_age / p_sex +
  plot_annotation(title = "Demographic Data Reporting by State/Province")

# Save
ggsave("figs/CAM_1000/reporting_rates_map.jpg", p_combined,
       width = 10, height = 18, dpi = 300)

print("Saved to figs/CAM_1000/reporting_rates_map.jpg")
