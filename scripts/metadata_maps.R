# TO-DO: A lot of changes are necessary
# Need to make this compatible with larger map areas (Canada/US/Mexico)
# Need to standardize color schemes
# Should this merge with seq_effort_plots.R?
# Time dependent?

library(tidyverse, quietly = "T ")
library(usmap, quietly = "T")
library(ggplot2, quietly = "T")
library(dbplyr, quietly = "T")
library(duckdb)

scenario <- "USA"
fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)
df_meta <- tbl(con,"metadata")

total_seq_n <- df_meta |> summarize(n()) |> collect()
state_summary <- df_meta |>
  mutate(age_na = as.numeric(is.na(age_adj))) |>
  mutate(location_na = as.numeric(is.na(location))) |>
  group_by(division) |>
  summarize(perc_age_valid = 1 - sum(age_na)/n(),
            perc_location_exists = 1 - sum(location_na)/n()) |>
  rename(state = division) |>
  collect()

plot_usmap(data = state_summary,values='perc_age_valid') +
  labs(
    title = "Valid Ages by State",
    legend = "Percent Valid Age"
  ) +
  scale_fill_gradient(breaks = c(0.2,0.4,0.6,0.8,1),
                       name="Percent Valid Age",
                       low = "white",
                       high = "orange"
  ) +
theme(legend.position=c(0.9,0.05))
ggsave("figs/metadata/map_state_perc_valid_age.jpg",dpi=300,width = 10.5, height = 5, units = "in")


plot_usmap(data = state_summary,values='perc_location_exists') +
  labs(
    title = "Location Data by State",
    legend = "Percent Location Included"
  ) +
  scale_fill_gradient(breaks = c(0.2,0.4,0.6,0.8,1),
                      name="Percent Location Included",
                      low = "white",
                      high = "orange"
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/metadata/map_state_perc_location_exists.jpg",dpi=300,width = 10.5, height = 5, units = "in")
