library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(argparse)
source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
SCENARIO <- args$scenario
fig_path <- paste0("figs/",SCENARIO,"/time/")

# Reference date for fold change calculations
REF_DATE <- as.Date("2023-01-01")

state_rr_snap <- read_tsv(paste0("results/",
                                 SCENARIO,
                                 "/time_state/df_state_rr_snap.tsv"))
state_rr_series <- read_tsv(paste0("results/",
                                   SCENARIO,
                                   "/time_state/df_state_rr_series.tsv"))
state_pair_types <- state_rr_series %>%
  select(x, y, pair_type) %>%
  distinct()
state_rr_snap <- state_rr_snap %>%
  left_join(state_pair_types, by = c("x", "y"))

df_regions <- read_csv("data/us_states_regions.csv") %>%
  select(state, bea_reg)

df_dist <- read_tsv("data/nb_dist_states.tsv")

PAIR_LIST <- tribble(
  ~x,                    ~y,
  "Mexico",              "Texas",
  "Mexico",              "California",
  "Mexico",              "Arizona",
  "British Columbia",    "Washington",
  "Ontario",             "Michigan",
  "Ontario",             "New York",
  "Ontario",             "Quebec",
  "Ontario",             "British Columbia",
  "Ontario",             "Alberta",
  "Quebec",              "New York",
  "Quebec",              "Vermont",
  "California",          "New York",
  "California",          "Arizona",
  "California",          "Texas",
  "California",          "Washington",
  "Michigan",            "Illinois",
  "Alberta",             "Montana",
  "Manitoba",            "North Dakota",
  "Manitoba",            "Alberta",
  "Vermont",             "New York" 
)

#Make quarter labels
format_quarters <- function(date){
  paste0(lubridate::year(date)," Q",lubridate::quarter(date))
}

ref_nRR <- state_rr_snap %>%
  filter(date == REF_DATE) %>%
  select(x,y,nRR) %>%
  rename(nRR_ref = nRR)


sub_RR_snap <- state_rr_snap %>%
  left_join(ref_nRR,by=c("x","y")) %>%
  inner_join(PAIR_LIST, by=c("x","y")) %>%
  mutate(pair_name = paste0(x,"-",y)) %>%
  mutate(is_intl = pair_type == "International") %>%
  # Create separate factor levels for each facet with alphabetical order (reversed for plotting)
  group_by(is_intl) %>%
  arrange(pair_name, .by_group = TRUE) %>%
  mutate(pair_name = factor(pair_name, levels = rev(unique(pair_name)))) %>%
  # Add fold change calculation
  group_by(pair_name) %>%
  arrange(date) %>%
  mutate(nRR_fold = nRR / nRR_ref) %>%
  ungroup()

# Create common scale for both plots
common_fill_scale <- scale_fill_distiller(
  palette = "RdYlBu",
  limits = c(-1, 1),
  oob = scales::squish,
  name = "Fold Change",
  labels = function(x) format(10^x, digits = 1)
)

# Create separate plots for each facet
p1 <- ggplot(sub_RR_snap %>% filter(!is_intl),
       aes(x=date, y=pair_name, fill = log10(nRR_fold))) +
  geom_tile(color="#555") +
  geom_vline(xintercept = as.numeric(REF_DATE),
             color = "black", linewidth = 1.5, linetype = "dotted") +
  scale_x_date(
    date_breaks = "3 months",
    labels = function(x) format_quarters(x),
    expand = expansion(0.01)
  ) +
  labs(x = "Date", y = "Domestic Pairs") +
  common_fill_scale +  # Use common scale
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    panel.grid.major = element_blank()
  )

p2 <- ggplot(sub_RR_snap %>% filter(is_intl),
       aes(x=date, y=pair_name, fill=log10(nRR_fold))) +
  geom_tile(color="#555") +
  geom_vline(xintercept = as.numeric(REF_DATE),
             color = "black", linewidth = 1.5, linetype = "dotted") +
  scale_x_date(
    date_breaks = "3 months",
    labels = function(x) format_quarters(x),
    expand = expansion(0.01)
  ) +
  labs(x = "Date", y = "International Pairs") +
  common_fill_scale +  # Use common scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8),
    legend.position = "none",
    panel.grid.major = element_blank()
  ) 
# Combine plots using patchwork
sub_pair_plot <- p1 / p2 +
  plot_annotation(
    title = "Normalized RR Over Time by State/Province Pair",
    theme = theme_minimal())

ggsave(filename = paste0(fig_path,"state_pair_nRR_heatmap.png"),plot = sub_pair_plot,
       width=8,
       height=5,
       dpi=192)


full_pair_snap <- state_rr_snap %>%
  filter(x != y) %>%
  left_join(ref_nRR, by = c("x", "y")) %>%
  mutate(pair_type = factor(pair_type,
                            levels = c("Different States, Same Region",
                                       "Different Regions, Same Country",
                                       "International"))) %>%
  mutate(pair_name = paste0(x,"-",y)) %>%
  arrange(pair_type, pair_name) %>%
  mutate(pair_name = factor(pair_name, levels = unique(pair_name))) %>%
  mutate(nRR_fold = nRR / nRR_ref) %>%
  ungroup()

# Create all pairs plot with faceting by pair_type
p_all_fold <- ggplot(full_pair_snap,
       aes(x=date, y=pair_name, fill = log10(nRR_fold))) +
  geom_tile(color="#555", linewidth=0) +
  scale_x_date(
    date_breaks = "3 months",
    labels = function(x) format_quarters(x),
    expand = expansion(0.01)
  ) +
  facet_grid(rows = vars(pair_type), scales = "free_y", space = "free_y") +
  common_fill_scale +
  labs(x = "Date", y = "State/Province Pairs") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8),
    axis.text.y = element_blank(),
    legend.position = "right",
    panel.grid.major = element_blank(),
    strip.text.y = element_text(angle = 270, hjust = 0.5)
  )

ggsave(paste0(fig_path,"all_pairs_fold_heatmap.png"),
       plot = p_all_fold,
       width=8,
       height=20,
       dpi=192)

full_pair_snap_dist <- full_pair_snap %>%
  left_join(df_dist,
     by=join_by("x" == "state_x","y" == "state_y")
  ) %>%
  mutate(pair_type = case_when(
    pair_type == "International" & nb_dist == 1 ~ "International, Bordering State/Provinces",
    pair_type == "International" & nb_dist != 1 ~ "International, Non-Bordering State/Provinces",
    pair_type == "International" & is.na(nb_dist) ~ "International, Non-Bordering State/Provinces",
    TRUE ~ as.character(pair_type)
))

gg_boxplot_fold <- ggplot(full_pair_snap_dist,
  aes(x = date,
    y = nRR_fold,
    color = pair_type,
    group=interaction(date,pair_type))
  )+
  geom_boxplot(outliers=FALSE) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "bottom") 

unique_series <- state_rr_series %>% 
  mutate(duplicate = x >= y) %>% 
  filter (duplicate != TRUE) %>%
  mutate(pair_name = paste0(x,"-",y)) %>%
  arrange(pair_name) %>%
  mutate(pair_type = factor(pair_type,
                            levels = c("Same State/Province",
                                       "Different States, Same Region",
                                       "Different Regions, Same Country",
                                       "International")))

ggplot(unique_series,
       aes(x = date, y = nRR, group = date, colour = pair_type)) +
  geom_boxplot(alpha = 0.2) +
  coord_cartesian(ylim=c(0,0.5)) +
  facet_wrap(~ pair_type, ncol = 1) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(paste0(fig_path,"rr_series.png"),
       width=8,
       height=8,
       dpi=192)

# Create poster version with all categories in one plot (2:1 aspect ratio, no outliers)
# Use rr_snap instead of rr_series to reduce clutter
unique_snap_filtered <- state_rr_snap %>%
  filter(pair_type != "Same State/Province", x != y) %>%
  mutate(pair_type = factor(pair_type,
                            levels = c("Different States, Same Region",
                                      "Different Regions, Same Country",
                                      "International")))

rr_series_poster <- ggplot(unique_snap_filtered,
       aes(x = date, y = nRR, group = interaction(date, pair_type), color = pair_type, fill = pair_type)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, linewidth = 0.4, position = position_dodge(width = 80)) +
  coord_cartesian(ylim = c(0, 0.4)) +
  scale_color_manual(
    values = c(
      "Different States, Same Region" = "lightcoral",
      "Different Regions, Same Country" = "#66C2A5",
      "International" = "cornflowerblue"
    ),
    name = "Pair Type"
  ) +
  scale_fill_manual(
    values = c(
      "Different States, Same Region" = "lightcoral",
      "Different Regions, Same Country" = "#66C2A5",
      "International" = "cornflowerblue"
    ),
    name = "Pair Type"
  ) +
  scale_x_date(date_breaks = "3 months",
               labels = function(x) format_quarters(x),
               expand = expansion(mult = 0.02)) +
  labs(x = "Date", y = "nRR", title = "Normalized Identical Sequence RR Over Time") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave(paste0(fig_path,"rr_series_poster.png"),
       plot = rr_series_poster,
       width = 10,
       height = 5,
       dpi = 192)

# ============================================================================
# Correlation time series poster (using snap data)
# ============================================================================

# Load travel data and distance data, join with snap data
df_travel <- read_tsv("data/travel_vars.tsv", show_col_types = FALSE)
df_cbsa_dist <- read_csv("data/state_pair_min_cbsa_distance_km.csv", show_col_types = FALSE) %>%
  rename(min_cbsa_dist = min_cbsa_distance_km)

# Function to normalize travel RR (same method as for nRR)
normalized_travel_rr <- function(df_rr){
  # Filter to the diagonal (RR(x,x)) entries
  rr_diag <- df_rr %>%
    filter(x == y) %>%
    select(date, state = x, RR_diag = RR) %>%
    mutate(RR_diag = ifelse(is.na(RR_diag), NA, RR_diag))

  # Join to get RR(x,x) and RR(y,y) for each row
  df_rr <- df_rr %>%
    left_join(rr_diag, by = c("date", "x" = "state")) %>%
    rename(RR_xx = RR_diag) %>%
    left_join(rr_diag, by = c("date", "y" = "state")) %>%
    rename(RR_yy = RR_diag) %>%
    mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
  return(df_rr)
}

# Join travel data and distance data, then normalize
state_rr_snap_travel <- state_rr_snap %>%
  left_join(df_travel, by = c("x", "y")) %>%
  left_join(df_cbsa_dist, by = c("x", "y"))

# Normalize RR_trips
if("RR_trips" %in% names(state_rr_snap_travel)) {
  df_trips <- state_rr_snap_travel %>%
    select(date, x, y, RR = RR_trips) %>%
    filter(!is.na(RR)) %>%
    normalized_travel_rr() %>%
    select(date, x, y, nRR_trips = nRR)
  state_rr_snap_travel <- left_join(state_rr_snap_travel, df_trips, by = c("date", "x", "y"))
}

# Normalize RR_move
if("RR_move" %in% names(state_rr_snap_travel)) {
  df_move <- state_rr_snap_travel %>%
    select(date, x, y, RR = RR_move) %>%
    filter(!is.na(RR)) %>%
    normalized_travel_rr() %>%
    select(date, x, y, nRR_move = nRR)
  state_rr_snap_travel <- left_join(state_rr_snap_travel, df_move, by = c("date", "x", "y"))
}

# Function to calculate correlations for snap data
calculate_snap_correlations <- function(data, var_list, var_labels) {
  results_list <- vector("list", length(var_list))

  for (i in seq_along(var_list)) {
    var <- var_list[i]
    var_data <- data %>% filter(x != y)

    dates <- unique(var_data$date) %>% sort()
    corr_results <- data.frame(
      date = dates,
      correlation = numeric(length(dates)),
      variable = var_labels[i]
    )

    for (j in seq_along(dates)) {
      date_data <- var_data[var_data$date == dates[j], ]
      valid_idx <- !is.na(date_data[[var]]) & !is.na(date_data$nRR)

      if (sum(valid_idx) >= 3) {
        corr_results$correlation[j] <- cor(
          date_data[[var]][valid_idx],
          date_data$nRR[valid_idx],
          method = "spearman"
        )
      } else {
        corr_results$correlation[j] <- NA
      }
    }

    results_list[[i]] <- corr_results
  }

  return(do.call(rbind, results_list))
}

# Calculate correlations for snap data
corr_snap <- calculate_snap_correlations(
  data = state_rr_snap_travel,
  var_list = c("min_cbsa_dist", "nRR_trips", "nRR_move"),
  var_labels = c("Geographic Distance", "Normalized Travel", "Normalized Mobility")
)

# Create correlation poster
correlation_poster <- ggplot(corr_snap, aes(x = date, y = correlation, color = variable)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c(
      "Geographic Distance" = "lightcoral",
      "Normalized Travel" = "#66C2A5",
      "Normalized Mobility" = "cornflowerblue"
    ),
    name = "Variable"
  ) +
  scale_x_date(date_breaks = "3 months",
               labels = format_quarters,
               expand = expansion(mult = 0.02)) +
  labs(
    title = "Spearman Correlation with Normalized RR Over Time",
    x = "Date",
    y = "Spearman Correlation"
  ) +
  ylim(-1, 1) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3)
  )

ggsave(paste0(fig_path,"correlations_poster.png"),
       plot = correlation_poster,
       width = 10,
       height = 5,
       dpi = 192)

ts_nRR_df <- unique_series %>%
  mutate(
    date_num  = as.numeric(date),
    pair_type = factor(pair_type),
    pair_name = factor(pair_name)
  ) %>%
  filter(!is.na(nRR), !is.na(date_num), !is.na(pair_type), !is.na(pair_name))


#ts_nRR_model <- mgcv::gam(
#  nRR ~ s(date_num,bs="cs",k=10) + pair_type + s(pair_name,bs="re"),
#  data=ts_nRR_df
#)

#summary(ts_nRR_model)
#mgcv::gam.check(ts_nRR_model)
#AIC(ts_nRR_model)

#Function that makes a 2x2 set of plots looking at self-enrichment vs pair-enrichment
#s1 and s2 are state names to be compared
state_pair_plot <- function(s1,s2){
  plot <- ggplot(data=(state_rr_series |> 
                         filter(x %in% c(s1,s2) & y %in% c(s1,s2))),
                 aes(x=date,y=RR,color=x)) +
    geom_line() +
    geom_point(aes(fill=x)) +
    geom_line(y=0,lty=2,size=1,color="black") +
    facet_grid(rows=vars(x),cols=vars(y)) + 
    scale_x_date(date_breaks = "6 months") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size= 11)) +
    theme(axis.text.y = element_text(size=11)) +
    labs(x = "Calendar Date", y = "RR",
         title=paste0("Monthly Enrichment - ",s1," & ",s2)) +
    theme(axis.title.x = element_text(size=14)) +
    theme(axis.title.y = element_text(size=14)) +
    theme(legend.position = "none")
  return(plot)
}

state_pair_plot("New York", "California")

state_list <- list("California",
                   "Washington",
                   "New York",
                   "New Jersey",
                   "Illinois",
                   "Wisconsin",
                   "Mexico",
                   "Ontario")

state_nrr_timeline <- function(state_name){
    nrr_series_sub <- state_rr_series %>%
        filter(x == state_name & y != state_name) %>%
        rename(paired_state = y) %>%
        select(date, paired_state, nRR) %>%
        left_join(df_regions, by = c("paired_state" = "state"))
        
    plot <- ggplot(nrr_series_sub,
                   aes(x = date, y = nRR, 
                       group = paired_state, 
                       color = bea_reg)) +
        geom_line(alpha = 0.5, linewidth = 0.3) +
        geom_smooth(aes(group=bea_reg), size=1) +
        scale_x_date(name = "Date",date_breaks = "6 months") +
        scale_y_continuous(name = "nRR") +
        ggtitle(paste0("nRR Time Series - ",state_name)) +
        region_color_scale() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    
    fn_plot <- paste0(fig_path,
                      "/individual_states/nRR_timeline_",
                      state_name,
                      ".png")
    ggsave(plot = plot,
        filename=fn_plot,
        width=8,
        height=5,
        dpi=192,
        create.dir = TRUE)
    return(plot)
}

lapply(state_list,state_nrr_timeline)
