# File: airline_rr_by_year.R
# Author(s): Amin Bemanian
# Date: 5/19/26
# Description: Calculate T100 airline RR by year and correlate against 2019 baseline

suppressPackageStartupMessages({
    library(tidyverse)
})

EXCL_STATES <- c("PR", "TT", "VI")
DATA_DIR    <- "data/T100_Dataset/T100D_Market_Year"
YEARS       <- 2019:2025
REF_YEAR    <- 2019

abbrev_to_state <- function(abbrev) {
    lookup <- c(state.abb, "DC")
    names  <- c(state.name, "District of Columbia")
    names[match(toupper(abbrev), lookup)]
}

calc_rr_t100 <- function(path) {
    raw <- read_csv(path, show_col_types = FALSE)

    agg <- raw %>%
        filter(!(ORIGIN_STATE_ABR %in% EXCL_STATES)) %>%
        filter(!(DEST_STATE_ABR   %in% EXCL_STATES)) %>%
        select(PASSENGERS, ORIGIN_STATE_ABR, DEST_STATE_ABR) %>%
        group_by(ORIGIN_STATE_ABR, DEST_STATE_ABR) %>%
        summarise(pass = sum(PASSENGERS), .groups = "drop") %>%
        rename(ori = ORIGIN_STATE_ABR, dest = DEST_STATE_ABR)

    uniq_states <- union(agg$ori, agg$dest)
    all_grid <- expand_grid(ori = uniq_states, dest = uniq_states, plus_one = 1)

    agg <- right_join(agg, all_grid, by = join_by(ori, dest)) %>%
        mutate(pass = ifelse(is.na(pass), 0, pass) + plus_one) %>%
        select(-plus_one)

    mirror <- agg %>% rename(ori_mirror = dest, dest_mirror = ori, pass_mirror = pass)

    left_join(agg, mirror, by = join_by("ori" == "ori_mirror", "dest" == "dest_mirror")) %>%
        mutate(
            x = abbrev_to_state(ori),
            y = abbrev_to_state(dest),
            pass_xy    = pass + pass_mirror,
            pass_total = sum(pass_xy)
        ) %>%
        group_by(x) %>% mutate(pass_x = sum(pass_xy)) %>% ungroup() %>%
        group_by(y) %>% mutate(pass_y = sum(pass_xy)) %>% ungroup() %>%
        mutate(RR_air = pass_xy * pass_total / pass_x / pass_y) %>%
        select(x, y, RR_air)
}

# Calculate RR for each year
message("Calculating RR by year...")
rr_by_year <- map(YEARS, function(yr) {
    path <- file.path(DATA_DIR, sprintf("T_T100D_MARKET_US_CARRIER_%d.csv", yr))
    message("  Processing ", yr)
    calc_rr_t100(path) %>% mutate(year = yr)
}) %>% bind_rows()

# Pivot wide for correlation analysis
rr_wide <- rr_by_year %>%
    filter(!is.na(x), !is.na(y)) %>%
    pivot_wider(names_from = year, values_from = RR_air, names_prefix = "RR_")

ref_col <- paste0("RR_", REF_YEAR)

# Pearson correlations against reference year (log10 scale, exclude zeros)
message("\nPearson correlations vs ", REF_YEAR, " (log10 scale):")
cor_results <- map_dfr(YEARS, function(yr) {
    col <- paste0("RR_", yr)
    df  <- rr_wide %>%
        filter(.data[[ref_col]] > 0, .data[[col]] > 0) %>%
        filter(!is.na(.data[[ref_col]]), !is.na(.data[[col]]))
    r <- cor(log10(df[[ref_col]]), log10(df[[col]]),
             method = "pearson", use = "complete.obs")
    n <- nrow(df)
    tibble(year = yr, r = round(r, 4), n_pairs = n)
})

print(cor_results)

# Save wide RR table
write_tsv(rr_wide, "data/rr_air_t100_by_year.tsv", na = "NA")
message("\nWide RR table saved to data/rr_air_t100_by_year.tsv")

# Scatter plot: each year vs 2019
p_scatter <- rr_by_year %>%
    filter(!is.na(x), !is.na(y)) %>%
    left_join(rr_by_year %>% filter(year == REF_YEAR) %>%
                  select(x, y, RR_ref = RR_air),
              by = c("x", "y")) %>%
    filter(year != REF_YEAR, RR_air > 0, RR_ref > 0) %>%
    left_join(cor_results %>% select(year, r), by = "year") %>%
    mutate(year_label = paste0(year, " (r = ", r, ")")) %>%
    ggplot(aes(x = RR_ref, y = RR_air)) +
    geom_point(alpha = 0.08, size = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linewidth = 0.8) +
    scale_x_continuous(transform = "log10",
                       name = paste0("T100 RR (", REF_YEAR, ")"),
                       breaks = 10^(-4:3),
                       labels = scales::label_log(base = 10),
                       limits = c(10^-4, 10^3)) +
    scale_y_continuous(transform = "log10",
                       name = "T100 RR (year)",
                       breaks = 10^(-4:3),
                       labels = scales::label_log(base = 10),
                       limits = c(10^-4, 10^3)) +
    facet_wrap(vars(year_label), ncol = 3) +
    theme_classic(base_size = 16) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_blank(),
          strip.text = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 13))

dir.create("figs/", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/rr_t100_by_year_vs_2019.jpg", plot = p_scatter,
       dpi = 300, width = 12, height = 8)
message("Scatter plot saved to figs/rr_t100_by_year_vs_2019.jpg")

# Monthly passengers from yearly files
pax_monthly <- map_dfr(YEARS, function(yr) {
    path <- file.path(DATA_DIR, sprintf("T_T100D_MARKET_US_CARRIER_%d.csv", yr))
    raw  <- read_csv(path, show_col_types = FALSE)
    raw %>%
        group_by(MONTH) %>%
        summarise(total_pax = sum(PASSENGERS, na.rm = TRUE), .groups = "drop") %>%
        mutate(date = as.Date(paste0(yr, "-", sprintf("%02d", MONTH), "-15")))
})

# Align correlation results to July 1
cor_plot_data <- cor_results %>%
    mutate(date = as.Date(paste0(year, "-07-01")))

# Scale factor for dual axis
pax_scale <- max(cor_plot_data$r) / max(pax_monthly$total_pax)

p_rho <- ggplot() +
    geom_bar(data = pax_monthly,
             aes(x = date, y = total_pax * pax_scale),
             stat = "identity", fill = "steelblue", alpha = 0.35, width = 25) +
    geom_line(data = cor_plot_data,
              aes(x = date, y = r), linewidth = 1, color = "black") +
    geom_point(data = cor_plot_data,
               aes(x = date, y = r), size = 4, shape = 21,
               fill = "firebrick", color = "black", stroke = 1.5) +
    geom_vline(xintercept = as.Date(paste0(REF_YEAR, "-07-01")),
               linetype = "dashed", color = "grey50") +
    scale_x_date(name = "", date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(
        name = "Pearson r vs 2019 (log10 scale)",
        limits = c(0, 1),
        sec.axis = sec_axis(~ . / pax_scale / 1e6,
                            name = "Monthly passengers (millions)")
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.y.right = element_text(color = "steelblue"),
          axis.text.y.right  = element_text(color = "steelblue"))

ggsave("figs/rr_t100_yearly_correlation.jpg", plot = p_rho,
       dpi = 300, width = 8, height = 4)
message("Correlation plot saved to figs/rr_t100_yearly_correlation.jpg")
