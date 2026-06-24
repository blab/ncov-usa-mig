#File: state_dist_plot.R
#Author(s): Amin Bemanian
#Date: 7/29/24
#Description: Makes a distance plots from state RR matrix
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggsignif)
library(patchwork)
library(scales)
library(readr)
library(cowplot)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--air_source', type = 'character', default = 'db1b',
                      choices = c('db1b', 't100'),
                      help = 'Air travel data source: db1b (full itinerary) or t100 (individual legs)')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
air_source <- args$air_source

####TEMP DELETE ONCE YOU WANT TO MAKE THE FLOW WORK
scenario <- "CAM_1000"

fig_path     <- paste0("figs/", scenario, "/dist/")
air_fig_path <- paste0("figs/", scenario, "/dist/air_", air_source, "/")
dir.create(fig_path,     recursive = TRUE, showWarnings = FALSE)
dir.create(air_fig_path, recursive = TRUE, showWarnings = FALSE)

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
state_rr <- fread(fn_rr)  %>%
  mutate(simple_adj = case_when(
    nb_dist == 0 ~ "Within state",
    nb_dist == 1 ~ "Adjacent state",
    TRUE ~ "Non-adjacent state"
  ))

state_rr$simple_adj<-factor(state_rr$simple_adj,levels=c("Non-adjacent state","Adjacent state","Within state"))

# Load travel data and select the appropriate air source column as RR_air
df_travel <- fread("data/travel_vars.tsv") %>%
  rename(RR_air = paste0("RR_air_", air_source))
state_rr <- left_join(state_rr, df_travel, by = c("x", "y"))

Y_AXIS_ZERO <- min(state_rr$RR) * 1E-2 #Arbitrarily lower value for log-transform

df_adj_boxplot <- state_rr %>% 
  group_by(simple_adj) %>% 
  summarise(y05 = quantile(RR, 0.05),
            y25 = quantile(RR, 0.25),
            y50 = quantile(RR, 0.5),
            y75 = quantile(RR, 0.75),
            y95 = quantile(RR, 0.95)) %>% 
  mutate(y05_plot = ifelse(y05 == 0., Y_AXIS_ZERO, y05), #Make sure no zero values
         y25_plot = ifelse(y25 == 0., Y_AXIS_ZERO, y25),
         y50_plot = ifelse(y50 == 0., Y_AXIS_ZERO, y50),
         y75_plot = ifelse(y75 == 0., Y_AXIS_ZERO, y75),
         y95_plot = ifelse(y95 == 0., Y_AXIS_ZERO, y95))

fn_state_adj_plot <- paste0(fig_path, "state_adj_plot.jpg")
state_adj_plot <- state_rr%>%
  ggplot(aes(x=simple_adj)) +
  geom_jitter(aes(y=RR),alpha=0.1,fill='black',size=1,width=0.3, height = 0) +
  geom_signif(aes(y=RR),
              comparisons = list(c("Within state","Adjacent state"), c("Adjacent state","Non-adjacent state")),
              map_signif_level = function(p) {
                ifelse(p < 0.001, "< 0.001", formatC(p, digits = 3, format = "f"))
              },
              tip_length = 0.01,
              textsize = 4,
              size = 1.5,
              test = "wilcox.test",
              y_position = log10(c(100, 10)), 
              color = 'firebrick') +
  geom_boxplot(data = df_adj_boxplot,
               aes(ymin = y05_plot, lower = y25_plot, middle = y50_plot,
                   upper = y75_plot, ymax = y95_plot),
               stat = "identity", fill = NA, width = 0.7) +
  scale_x_discrete(name= '') +
  scale_y_continuous(transform ='log10',
                     name="seqRR",
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2, 1E3),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2}),expression(10^{3})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-1.5),10^(3))) +
  coord_flip() +
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  ggtitle(scenario)  
state_adj_plot

ggsave(fn_state_adj_plot,
       plot=state_adj_plot,
       device = "jpeg",
       dpi = 300,
       width = 8,
       height = 5
)

rho_euclid <- cor(state_rr$euclid_dist, log(state_rr$RR), method = "pearson", use = "complete.obs")
rho_cbsa <- cor(state_rr$min_cbsa_dist, log(state_rr$RR), method = "pearson", use = "complete.obs")


fn_state_euclid_dist_plot <- paste0(fig_path, "state_euclid_dist_plot.jpg")
state_euclid_dist_plot <- state_rr %>%
  filter(x != y) %>%
  ggplot() +
  geom_point(aes(x=euclid_dist,y=RR),alpha=0.1,size=1,fill='black') +
  geom_smooth(aes(x=euclid_dist,y=RR),color='firebrick',fill='firebrick',method='loess') +
  scale_x_continuous(name="State Centroid Distance (km)",
                     breaks=c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),
                     limits=c(0,2000)) +
  scale_y_continuous(transform ='log',
                     name="seqRR",
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.15),10^(1.5))) +
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(scenario) +
  labs(subtitle = paste0("r = ",round(rho_euclid,2)))

ggsave(fn_state_euclid_dist_plot,
       plot=state_euclid_dist_plot,
       device = "jpeg",
       dpi = 300,
       width = 6,
       height = 4
)

fn_state_cbsa_dist_plot <- paste0(fig_path, "state_cbsa_dist_plot.jpg")
state_cbsa_dist_plot <- state_rr %>%
  filter(x != y) %>%
  ggplot() +
  geom_point(aes(x=min_cbsa_dist,y=RR),alpha=0.1,size=1,fill='black') +
  geom_smooth(aes(x=min_cbsa_dist,y=RR),color='firebrick',fill='firebrick',method='loess') +
  scale_x_continuous(name="State Nearest CBSA Distance (km)",
                     breaks=c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),
                     limits=c(0,2000)) +
  scale_y_continuous(transform ='log',
                     name="seqRR",
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.15),10^(1.5))) +
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(scenario) +
  labs(subtitle = paste0("r = ",round(rho_cbsa,2)))

ggsave(fn_state_cbsa_dist_plot,
       plot=state_cbsa_dist_plot,
       device = "jpeg",
       dpi = 300,
       width = 6,
       height = 4
)

fn_euclid_cbsa_dist_plot <- paste0(fig_path, "state_euclid_vs_cbsa.jpg")
state_euclid_cbsa_dist_plot <- state_rr %>%
  filter(x != y) %>%
  ggplot() +
  geom_point(aes(x=euclid_dist,y=RR),alpha = 0.1,color="firebrick") + 
  geom_smooth(aes(x=euclid_dist,y=RR),alpha=0.5,color = "firebrick") + 
  geom_point(aes(x=min_cbsa_dist,y=RR),alpha=0.1,color="navyblue") + 
  geom_smooth(aes(x=min_cbsa_dist,y=RR),alpha=0.5,color = "navyblue") + 
  scale_x_continuous(name="Distance (km)",
                     breaks=c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),
                     limits=c(0,3000)) +
  scale_y_continuous(transform ='log',
                     name="seqRR",
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.5),10^(1.5))) +
  geom_hline(yintercept = 1) + #Reference point
  theme_bw()

fn_state_euclid_logdist_plot <- paste0(fig_path, "state_euclid_logdist_plot.jpg")
state_euclid_logdist_plot <- state_rr %>%
  ggplot() +
  geom_point(aes(x=log10(euclid_dist+1),y=RR),alpha=0.1,size=1,fill='black') +
  geom_smooth(aes(x=log10(euclid_dist+1),y=RR),color='firebrick',fill='firebrick',method='loess') +
  scale_x_continuous(name=expression(log["10"]("Centroid Distance")),
                     breaks = c(0,1,2,3,4,5),
                     labels = c(0,expression(10^{1}),expression(10^{2}),expression(10^{3}),expression(10^{4}),expression(10^{5})),
                     limits=c(0,4)) +
  scale_y_continuous(transform ='log',
                     name="seqRR",
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.15),10^(1.5))) +
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(scenario) +
  labs(subtitle = paste0("r = ",round(rho_euclid,2)))

ggsave(fn_state_euclid_logdist_plot,
       plot=state_euclid_logdist_plot,
       device = "jpeg",
       dpi = 300,
       width = 6,
       height = 6
)

fn_state_nb_dist_plot <- paste0(fig_path, "state_nb_dist_plot.jpg")
state_nb_dist_plot <- state_rr %>%
  ggplot() +
  geom_boxplot(aes(x = nb_dist, y = RR, group = nb_dist),
             fill = 'black', linewidth = 0.5) +
  geom_smooth(aes(x=nb_dist,y=RR),color='firebrick',fill='firebrick',method='loess') +
  scale_x_continuous(name = "Neighbor Order (Queen-Adjacency)",
                     breaks = 0:11,
                     labels = c("Self","1st","2nd","3rd","4th","5th","6th",
                     "7th","8th","9th","10th","11th"),
                     limits= c(0,4)) +
  scale_y_continuous(transform ='log',
                     name="seqRR",
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.15),10^(1.5))) +
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ggtitle(scenario)

ggsave(fn_state_nb_dist_plot,
       plot=state_nb_dist_plot,
       device = "jpeg",
       dpi = 300,
       width = 6,
       height = 4,
       create.dir = TRUE
)

# Create 2x2 patchwork of key distance/travel metrics
REGION_DATA <- fread("data/us_states_regions.csv")

# Add pair type: Same state / Intra-regional / Inter-regional
state_rr <- state_rr %>%
  left_join(REGION_DATA %>% select(state, bea_reg), by = c("x" = "state")) %>%
  rename(bea_reg_x = bea_reg) %>%
  left_join(REGION_DATA %>% select(state, bea_reg), by = c("y" = "state")) %>%
  rename(bea_reg_y = bea_reg) %>%
  mutate(pair_type = factor(
    case_when(
      x == y                             ~ "Same state",
      bea_reg_x == bea_reg_y             ~ "Intra-regional",
      TRUE                               ~ "Inter-regional"
    ),
    levels = c("Same state", "Intra-regional", "Inter-regional")
  ))

PAIR_COLORS <- c("Same state" = "cornflowerblue",
                 "Intra-regional" = "lightcoral",
                 "Inter-regional" = "lightgoldenrod")

# Calculate correlations. Exclude pairs with zero underlying movement so the
# floored/+correction points (esp. DB1B, where most pairs have no sampled
# itinerary) don't dominate the correlation.
air_count <- paste0("pass_xy_", air_source)

rho_cbsa <- state_rr %>%
  filter(x != y, RR > 0) %>%
  summarize(rho = cor(min_cbsa_dist, log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)
rho_trips <- state_rr %>%
  filter(x != y, trips_xy > 0, RR_trips > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_trips), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)
rho_move <- state_rr %>%
  filter(x != y, n_move_avg > 0, RR_move > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_move), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

# Calculate air travel correlation (zero-flight pairs excluded)
rho_air <- state_rr %>%
  filter(x != y, !is.na(RR_air), RR_air > 0, RR > 0, .data[[air_count]] > 2) %>%
  summarize(rho = cor(log(RR_air), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

POSTER_BASE  <- 12
POSTER_PT_SZ <- 0.6
POSTER_LW    <- 1.0
POSTER_SUBT  <- 11

poster_y_scale <- scale_y_continuous(
  transform = 'log',
  name = "seqRR",
  breaks = c(1E-1, 1, 1E1),
  labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
  limits = c(10^(-1.5), 10^(1.5))
)
poster_theme <- theme_classic(base_size = POSTER_BASE) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = POSTER_SUBT),
        legend.position = "none")

# 1. CBSA distance
p_cbsa <- state_rr %>%
  filter(x != y) %>%
  ggplot(aes(x = min_cbsa_dist, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = POSTER_PT_SZ) +
  geom_smooth(method = 'loess', linewidth = POSTER_LW, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = TRUE) +
  scale_x_continuous(name = "CBSA Distance (km)", limits = c(0, 3000)) +
  poster_y_scale +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  poster_theme +
  labs(subtitle = paste0("r = ", round(rho_cbsa, 2)))

# 2. Air Travel RR (DB1B)
p_air <- state_rr %>%
  filter(x != y, !is.na(RR_air), RR_air > 0, .data[[air_count]] > 2) %>%
  ggplot(aes(x = RR_air, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = POSTER_PT_SZ) +
  geom_smooth(method = 'loess', linewidth = POSTER_LW, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = paste0("Air Travel RR (", toupper(air_source), ")"),
                     transform = 'log10',
                     breaks = c(0.001, 0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-3}), expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  poster_y_scale +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  poster_theme +
  labs(subtitle = paste0("r = ", round(rho_air, 2)))

# 3. Ground Travel RR (NHTS)
p_trips <- state_rr %>%
  filter(x != y, !is.na(RR_trips), trips_xy > 0) %>%
  ggplot(aes(x = RR_trips, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = POSTER_PT_SZ) +
  geom_smooth(method = 'loess', linewidth = POSTER_LW, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Ground Travel RR (NHTS)",
                     transform = 'log10',
                     breaks = c(0.0001, 0.01, 1, 100),
                     labels = c(expression(10^{-4}), expression(10^{-2}),
                               expression(10^{0}), expression(10^{2}))) +
  poster_y_scale +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  poster_theme +
  labs(subtitle = paste0("r = ", round(rho_trips, 2)))

# 4. Mobility RR (SafeGraph)
p_move <- state_rr %>%
  filter(x != y, !is.na(RR_move), n_move_avg > 0) %>%
  ggplot(aes(x = RR_move, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = POSTER_PT_SZ) +
  geom_smooth(method = 'loess', linewidth = POSTER_LW, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Mobility RR (SafeGraph)",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}), expression(10^{2}))) +
  poster_y_scale +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  poster_theme +
  labs(subtitle = paste0("r = ", round(rho_move, 2)))

shared_legend <- cowplot::get_legend(
  p_cbsa +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme(legend.position = "right", legend.text = element_text(size = POSTER_BASE - 1))
)
panels <- p_cbsa | p_air | p_trips | p_move
combined_dist_plot <- cowplot::plot_grid(panels, shared_legend, nrow = 1, rel_widths = c(1, 0.15))

fn_combined_dist <- paste0(fig_path, "state_distance_poster.jpg")
ggsave(fn_combined_dist,
       plot = combined_dist_plot,
       device = "jpeg",
       dpi = 300,
       width = 10,
       height = 2.5,
       units = "in",
       create.dir = TRUE
)
ggsave(paste0(fig_path, "state_distance_poster.svg"),
       plot = combined_dist_plot,
       device = "svg",
       width = 10,
       height = 2.5,
       units = "in"
)

# ============================================================================
# Air Travel RR by Flight Distance (3x1 vertical panel)
# ============================================================================

# Filter to different states only and where air travel data exists
state_rr_air <- state_rr %>%
  filter(x != y, !is.na(RR_air), !is.na(min_cbsa_dist))

# Stratify by flight distance bins (based on CBSA distance)
state_rr_air <- state_rr_air %>%
  mutate(flight_length = case_when(
    min_cbsa_dist < 300 ~ "Short",
    min_cbsa_dist > 1800 ~ "Long",
    TRUE ~ "Medium"
  )) %>%
  mutate(flight_length = factor(flight_length, levels = c("Short", "Medium", "Long")))

# Calculate correlations for each flight length bin
rho_short_air <- state_rr_air %>%
  filter(flight_length == "Short", RR_air > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_air), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

rho_medium_air <- state_rr_air %>%
  filter(flight_length == "Medium", RR_air > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_air), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

rho_long_air <- state_rr_air %>%
  filter(flight_length == "Long", RR_air > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_air), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

# Create plot for short flights
p_air_short <- state_rr_air %>%
  filter(flight_length == "Short") %>%
  ggplot(aes(x = RR_air, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = "seqRR",
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Short distance (<300 km)",
       subtitle = paste0("r = ", round(rho_short_air, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "none")

# Create plot for medium flights
p_air_medium <- state_rr_air %>%
  filter(flight_length == "Medium") %>%
  ggplot(aes(x = RR_air, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = "seqRR",
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Medium distance (300-1800 km)",
       subtitle = paste0("r = ", round(rho_medium_air, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "none")

# Create plot for long flights
p_air_long <- state_rr_air %>%
  filter(flight_length == "Long") %>%
  ggplot(aes(x = RR_air, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = "seqRR",
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Long distance (>1800 km)",
       subtitle = paste0("r = ", round(rho_long_air, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "none")

# Simple unstratified air travel correlation (reference)
rho_all_air <- state_rr_air %>%
  filter(RR_air > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_air), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

p_air_all <- state_rr_air %>%
  ggplot(aes(x = RR_air, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = "seqRR",
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "All distances",
       subtitle = paste0("r = ", round(rho_all_air, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

ggsave(paste0(air_fig_path, "state_air_all_dist.jpg"),
       plot = p_air_all, device = "jpeg", dpi = 300, width = 5, height = 4)

# Combine into 3x1 vertical grid
combined_air_plot <- p_air_short / p_air_medium / p_air_long

fn_combined_air <- paste0(air_fig_path, "state_air_distance_poster.jpg")
ggsave(fn_combined_air,
       plot = combined_air_plot,
       device = "jpeg",
       dpi = 300,
       width = 5,
       height = 8,
       create.dir = TRUE
)
ggsave(paste0(air_fig_path, "state_air_distance_poster.svg"),
       plot = combined_air_plot,
       device = "svg",
       width = 5,
       height = 8
)

message(paste("Air travel distance poster saved to:", fn_combined_air))

# ============================================================================
# Air Travel RR by Inter vs Intra-Regional (alternate stratification)
# ============================================================================

state_rr_air_reg <- state_rr_air %>%
  filter(!is.na(bea_reg_x), !is.na(bea_reg_y)) %>%
  mutate(flight_type = ifelse(bea_reg_x == bea_reg_y, "Intra-regional", "Inter-regional"),
         flight_type = factor(flight_type, levels = c("Intra-regional", "Inter-regional")))

rho_intra <- state_rr_air_reg %>%
  filter(flight_type == "Intra-regional", RR_air > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_air), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

rho_inter <- state_rr_air_reg %>%
  filter(flight_type == "Inter-regional", RR_air > 0, RR > 0) %>%
  summarize(rho = cor(log(RR_air), log(RR), method = "pearson", use = "complete.obs")) %>%
  pull(rho)

p_air_intra <- state_rr_air_reg %>%
  filter(flight_type == "Intra-regional") %>%
  ggplot(aes(x = RR_air, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = "seqRR",
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Intra-regional flights",
       subtitle = paste0("r = ", round(rho_intra, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

ggsave(paste0(air_fig_path, "state_air_intra_regional.jpg"),
       plot = p_air_intra, device = "jpeg", dpi = 300, width = 5, height = 4)

p_air_inter <- state_rr_air_reg %>%
  filter(flight_type == "Inter-regional") %>%
  ggplot(aes(x = RR_air, y = RR)) +
  geom_point(aes(color = pair_type), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = PAIR_COLORS, name = "", drop = FALSE) +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = "seqRR",
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Inter-regional flights",
       subtitle = paste0("r = ", round(rho_inter, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

ggsave(paste0(air_fig_path, "state_air_inter_regional.jpg"),
       plot = p_air_inter, device = "jpeg", dpi = 300, width = 5, height = 4)

# ============================================================================
# Conserved vs non-conserved inter-regional pairs: travel RR boxplots
# ============================================================================

get_conserved_pairs <- function(pct, scenario) {
  fread(paste0("results/", scenario, "/time_state/df_conserved_connections_", pct, "pct_5plus.tsv")) %>%
    mutate(edge_pair = paste(pmin(x, y), pmax(x, y), sep = "_")) %>%
    pull(edge_pair) %>%
    unique()
}

conserved_pairs_95 <- get_conserved_pairs(95, scenario)
conserved_pairs_90 <- get_conserved_pairs(90, scenario)
conserved_pairs_80 <- get_conserved_pairs(80, scenario)
conserved_pairs_70 <- get_conserved_pairs(70, scenario)

TIER_LEVELS  <- c(">=95th", "90-94th", "80-89th", "70-79th", "<70th")
TIER_COLORS  <- c("firebrick", "darkorange", "steelblue", "mediumpurple", "grey85")

assign_tier <- function(edge_pair) {
  factor(
    case_when(
      edge_pair %in% conserved_pairs_95 ~ ">=95th",
      edge_pair %in% conserved_pairs_90 ~ "90-94th",
      edge_pair %in% conserved_pairs_80 ~ "80-89th",
      edge_pair %in% conserved_pairs_70 ~ "70-79th",
      TRUE                               ~ "<70th"
    ),
    levels = TIER_LEVELS
  )
}

# sig_connections for quarterly conserved labels
sig_connections <- fread(paste0("results/", scenario, "/time_state/df_significant_connections_98_percentile.tsv")) %>%
  mutate(edge_pair = paste(pmin(x, y), pmax(x, y), sep = "_"),
         quarter = case_when(
           as.integer(format(as.Date(date), '%m')) %in% 1:3  ~ "Q1 (Jan-Mar)",
           as.integer(format(as.Date(date), '%m')) %in% 4:6  ~ "Q2 (Apr-Jun)",
           as.integer(format(as.Date(date), '%m')) %in% 7:9  ~ "Q3 (Jul-Sep)",
           TRUE                                               ~ "Q4 (Oct-Dec)"
         )) %>%
  filter(x < y, region_x != region_y)

conserved_by_quarter <- sig_connections %>%
  group_by(quarter, edge_pair) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 2) %>%
  select(quarter, edge_pair)

# Base inter-regional pairs with all travel metrics, deduplicated
inter_reg_base <- state_rr_air_reg %>%
  filter(flight_type == "Inter-regional", x < y) %>%
  mutate(edge_pair = paste(pmin(x, y), pmax(x, y), sep = "_"),
         conserved = assign_tier(edge_pair))

ADJACENT_COMPS <- list(c(">=95th", "90-94th"), c("90-94th", "80-89th"), c("80-89th", "70-79th"), c("70-79th", "<70th"))

make_conserved_boxplot <- function(data, y_var, y_label, title) {
  kw_p <- kruskal.test(reformulate("conserved", y_var),
                       data = data %>% filter(!is.na(.data[[y_var]]), .data[[y_var]] > 0))$p.value
  kw_label <- ifelse(kw_p < 0.001, "Kruskal-Wallis p < 0.001",
                     paste0("Kruskal-Wallis p = ", formatC(kw_p, digits = 3, format = "f")))
  data %>%
    filter(!is.na(.data[[y_var]]), .data[[y_var]] > 0) %>%
    ggplot(aes(x = conserved, y = .data[[y_var]], color = conserved)) +
    geom_boxplot(fill = NA, linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(alpha = 0.2, size = 1, width = 0.25, height = 0) +
    geom_signif(comparisons = ADJACENT_COMPS,
                map_signif_level = function(p) ifelse(p < 0.001, "***",
                                             ifelse(p < 0.01,  "**",
                                             ifelse(p < 0.05,  "*", "ns"))),
                tip_length = 0.01, textsize = 3.5, size = 0.5, step_increase = 0.08,
                test = "wilcox.test", color = "black") +
    scale_color_manual(values = TIER_COLORS, guide = "none") +
    scale_y_continuous(transform = 'log10', name = y_label,
                       breaks = c(0.01, 0.1, 1, 10),
                       labels = c(expression(10^{-2}), expression(10^{-1}),
                                 expression(10^{0}), expression(10^{1})),
                       expand = expansion(mult = c(0.05, 0.2))) +
    scale_x_discrete(name = "") +
    theme_classic(base_size = 11) +
    labs(title = title, subtitle = kw_label) +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
}

make_conserved_quarterly_boxplot <- function(base_data, y_var, y_label, title) {
  tidyr::crossing(
      edge_pair = base_data$edge_pair,
      quarter   = c("Q1 (Jan-Mar)", "Q2 (Apr-Jun)", "Q3 (Jul-Sep)", "Q4 (Oct-Dec)")
    ) %>%
    left_join(base_data %>% select(edge_pair, conserved, all_of(y_var)), by = "edge_pair") %>%
    mutate(quarter = factor(quarter, levels = c("Q1 (Jan-Mar)", "Q2 (Apr-Jun)",
                                                "Q3 (Jul-Sep)", "Q4 (Oct-Dec)"))) %>%
    filter(!is.na(.data[[y_var]]), .data[[y_var]] > 0) %>%
    ggplot(aes(x = conserved, y = .data[[y_var]], color = conserved)) +
    geom_boxplot(fill = NA, linewidth = 0.5, outlier.shape = NA) +
    geom_jitter(alpha = 0.2, size = 0.8, width = 0.25, height = 0) +
    geom_signif(comparisons = ADJACENT_COMPS,
                map_signif_level = function(p) ifelse(p < 0.001, "***",
                                             ifelse(p < 0.01,  "**",
                                             ifelse(p < 0.05,  "*", "ns"))),
                tip_length = 0.01, textsize = 3, size = 0.5, step_increase = 0.08,
                test = "wilcox.test", color = "black") +
    scale_color_manual(values = TIER_COLORS, guide = "none") +
    scale_y_continuous(transform = 'log10', name = y_label,
                       breaks = c(0.01, 0.1, 1, 10),
                       labels = c(expression(10^{-2}), expression(10^{-1}),
                                 expression(10^{0}), expression(10^{1})),
                       expand = expansion(mult = c(0.05, 0.2))) +
    scale_x_discrete(name = "") +
    facet_wrap(~quarter, nrow = 2) +
    theme_classic(base_size = 10) +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# Air travel RR (DB1B)
p_conserved_air <- make_conserved_boxplot(inter_reg_base, "RR_air", "Air Travel RR",
  "Conserved Percentile Pairs vs Air Travel RR (DB1B)")
ggsave(paste0(air_fig_path, "state_air_conserved_boxplot.jpg"),
       plot = p_conserved_air, device = "jpeg", dpi = 300, width = 5, height = 5)
ggsave(paste0(fig_path, "state_air_db1b_conserved_boxplot.jpg"),
       plot = p_conserved_air, device = "jpeg", dpi = 300, width = 5, height = 5)
ggsave(paste0(fig_path, "state_air_db1b_conserved_boxplot.svg"),
       plot = p_conserved_air, device = "svg", width = 7, height = 4)

p_conserved_air_quarterly <- make_conserved_quarterly_boxplot(inter_reg_base, "RR_air", "Air Travel RR",
  "Conserved Percentile Pairs vs Air Travel RR (DB1B) by Quarter")
ggsave(paste0(air_fig_path, "state_air_conserved_quarterly_boxplot.jpg"),
       plot = p_conserved_air_quarterly, device = "jpeg", dpi = 300, width = 7, height = 6)

# NHTS trips RR
p_conserved_trips <- make_conserved_boxplot(inter_reg_base, "RR_trips", "Travel RR (NHTS)",
  "Conserved Percentile Pairs vs NHTS Trips RR")
ggsave(paste0(fig_path, "state_nhts_conserved_boxplot.jpg"),
       plot = p_conserved_trips, device = "jpeg", dpi = 300, width = 5, height = 5)
ggsave(paste0(fig_path, "state_nhts_conserved_boxplot.svg"),
       plot = p_conserved_trips, device = "svg", width = 7, height = 4)

p_conserved_trips_quarterly <- make_conserved_quarterly_boxplot(inter_reg_base, "RR_trips", "Travel RR (NHTS)",
  "Conserved Percentile Pairs vs NHTS Trips RR by Quarter")
ggsave(paste0(fig_path, "state_nhts_conserved_quarterly_boxplot.jpg"),
       plot = p_conserved_trips_quarterly, device = "jpeg", dpi = 300, width = 7, height = 6)

# SafeGraph mobility RR
p_conserved_move <- make_conserved_boxplot(inter_reg_base, "RR_move", "Mobility RR (SafeGraph)",
  "Conserved Percentile Pairs vs SafeGraph Mobility RR")
ggsave(paste0(fig_path, "state_safegraph_conserved_boxplot.jpg"),
       plot = p_conserved_move, device = "jpeg", dpi = 300, width = 5, height = 5)
ggsave(paste0(fig_path, "state_safegraph_conserved_boxplot.svg"),
       plot = p_conserved_move, device = "svg", width = 7, height = 4)

p_conserved_move_quarterly <- make_conserved_quarterly_boxplot(inter_reg_base, "RR_move", "Mobility RR (SafeGraph)",
  "Conserved Percentile Pairs vs SafeGraph Mobility RR by Quarter")
ggsave(paste0(fig_path, "state_safegraph_conserved_quarterly_boxplot.jpg"),
       plot = p_conserved_move_quarterly, device = "jpeg", dpi = 300, width = 7, height = 6)
