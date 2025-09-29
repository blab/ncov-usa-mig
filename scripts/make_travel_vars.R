# File: make_travel_vars.R
# Author(s): Amin Bemanian
# Date: 9/25/25
# Description: Process travel and mobility covariates from different data sources
#            and combine them into a single dataset for analysis
# Arguments: None
# Output: data/travel_vars.tsv containing processed travel variables

options(warn = -1) # Suppress warnings
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
})

df_travel <- data.table::fread("data/interstate_travel_dot.csv") %>%
    rename(x = Origin) %>%
    rename(y = Destination) %>%
    rename(mode = "Pivot Field Names") %>%
    rename(trips_dir = "Pivot Field Values") %>% # Directional trips
    mutate(trips_dir = as.numeric(trips_dir)) # Switch from into to num since it will confuse the math later

STATE_LIST <- unique(df_travel$x)

# Convert from unidirectional to bidirectional trips
df_travel$trips_xy <- 0
for (i in STATE_LIST) {
    for (j in STATE_LIST) {
        if (i == j) {
            df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir
        } else {
            df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir +
                df_travel[x == j & y == i]$trips_dir
        }
    }
}

df_travel <- df_travel %>%
    group_by(x) %>%
    mutate(trips_x = sum(trips_xy)) %>%
    group_by(y) %>%
    mutate(trips_y = sum(trips_xy)) %>%
    ungroup() %>%
    mutate(
        trips_total = sum(trips_xy),
        same_state = (x == y), # Use to stratify for out of state travel only
        RR_trips = (trips_xy / trips_x) / (trips_y / trips_total)
    ) %>%
    select(x, y, same_state, RR_trips)

# Airline data already processed by calculate_airline_RR.R
df_travel_air <- read_csv("data/rr_air_db1b.csv")

df_move <- read_tsv("data/safegraph_states_adj_pullano.tsv") %>%
    rename(
        x = state_origin,
        y = state_destination,
        n_move_xy = n_movement_origin_destination,
    ) %>%
    mutate(n_move_all = sum(n_move_xy, na.rm = TRUE)) %>%
    group_by(x) %>%
    mutate(n_move_x = sum(n_move_xy, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(y) %>%
    mutate(n_move_y = sum(n_move_xy, na.rm = TRUE)) %>%
    ungroup()
df_move_rev <- df_move %>%
    select(x_rev = x, y_rev = y, n_move_xy) %>%
    rename(n_move_yx = n_move_xy)

df_move <- df_move %>%
    left_join(df_move_rev, by = c("x" = "y_rev", "y" = "x_rev")) %>%
    mutate(n_move_avg = 1 / 2 * (n_move_xy + n_move_yx)) %>%
    mutate(RR_move = n_move_avg * n_move_all / (n_move_x * n_move_y)) %>%
    select(x, y, RR_move)

df_travel <- df_travel %>%
    left_join(df_travel_air, by = c("x", "y")) %>%
    left_join(df_move, by = c("x", "y"))

write_tsv(df_travel, "data/travel_vars.tsv", na = "NA")
