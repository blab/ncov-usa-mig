# File: state_air_rr_map.R
# Author(s): Amin Bemanian
# Date: 5/20/26
# Description: Generate T100 and DB1B air travel RR maps for a target state
#              vs all other continental US states
# Arguments:
#   --state: Target state name (e.g. "New York")
#   --scenario: Scenario name (default: CAM_1000)
#   --rr_min: Lower bound of color scale (default: 0.5)
#   --rr_max: Upper bound of color scale (default: 2.0)

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(argparse)
    library(usmap)
    library(sf)
    library(patchwork)
})

collect_args <- function() {
    parser <- ArgumentParser()
    parser$add_argument("--state",    type = "character", required = TRUE,
                        help = "Target state name (e.g. 'New York')")
    parser$add_argument("--scenario", type = "character", default = "CAM_1000",
                        help = "Scenario name")
    parser$add_argument("--rr_min",   type = "double",    default = 0.5,
                        help = "Lower bound of RR color scale")
    parser$add_argument("--rr_max",   type = "double",    default = 2.0,
                        help = "Upper bound of RR color scale")
    parser$parse_args()
}

args     <- collect_args()
TARGET   <- args$state
scenario <- args$scenario
SCALE_LIM <- c(log10(args$rr_min), log10(args$rr_max))
EXCL     <- c("Alaska", "Hawaii", "Puerto Rico")

# Validate state name
valid_states <- c(state.name, "District of Columbia")
if (!(TARGET %in% valid_states)) {
    stop(sprintf("'%s' is not a recognized US state. Check spelling and capitalization.", TARGET))
}

df <- fread("data/travel_vars.tsv") %>%
    filter((x == TARGET | y == TARGET), x != y) %>%
    mutate(state = ifelse(x == TARGET, y, x)) %>%
    select(state, RR_air_t100, RR_air_db1b)

map_sf <- us_map("states") %>%
    filter(!(full %in% EXCL)) %>%
    left_join(df, by = c("full" = "state"))

make_map <- function(map_sf, rr_col, title, scale_lim) {
    ggplot(map_sf) +
        geom_sf(aes(fill = ifelse(full == TARGET, NA, log10(.data[[rr_col]]))),
                color = "white", linewidth = 0.3) +
        geom_sf(data = map_sf %>% filter(full == TARGET),
                fill = "grey40", color = "white", linewidth = 0.3) +
        scale_fill_gradient2(
            low      = "steelblue",
            mid      = "lightyellow",
            high     = "firebrick",
            midpoint = 0,
            na.value = "grey40",
            limits   = scale_lim,
            oob      = scales::squish,
            name     = "Air Travel RR",
            breaks   = log10(c(args$rr_min, 1, args$rr_max)),
            labels   = c(args$rr_min, "1", args$rr_max)
        ) +
        theme_void(base_size = 12) +
        labs(title = title) +
        theme(plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
              legend.position  = "bottom",
              legend.title     = element_text(size = 10),
              legend.key.width = unit(1.5, "cm"))
}

p_t100 <- make_map(map_sf, "RR_air_t100",
                   paste0("T100 Air Travel RR vs ", TARGET), SCALE_LIM)
p_db1b <- make_map(map_sf, "RR_air_db1b",
                   paste0("DB1B Air Travel RR vs ", TARGET), SCALE_LIM)

combined <- p_t100 / p_db1b + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

state_slug <- tolower(gsub(" ", "_", TARGET))
fn_out <- paste0("figs/", scenario, "/", state_slug, "_air_rr_map.jpg")
dir.create(dirname(fn_out), recursive = TRUE, showWarnings = FALSE)

ggsave(fn_out, plot = combined, dpi = 300, width = 8, height = 9)
message("Saved to: ", fn_out)
