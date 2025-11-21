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

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

####TEMP DELETE ONCE YOU WANT TO MAKE THE FLOW WORK
scenario <- "CAM_1000"

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
state_rr <- fread(fn_rr)  %>%
  mutate(simple_adj = case_when(
    nb_dist == 0 ~ "Within state",
    nb_dist == 1 ~ "Adjacent state",
    TRUE ~ "Non-adjacent state"
  ))

state_rr$simple_adj<-factor(state_rr$simple_adj,levels=c("Non-adjacent state","Adjacent state","Within state"))

# Load travel data
df_travel <- fread("data/travel_vars.tsv")
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

fn_state_adj_plot <- paste0("figs/",scenario,"/state_adj_plot.jpg") 
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
                     name=expression(RR["identical sequences"]),
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
       dpi = 200,
       width = 8,
       height = 5
)

#Note since spearman is rank-order correlation, log transformation does not matter
rho_euclid <- cor(state_rr$euclid_dist,state_rr$RR,method = "spearman",use = "complete.obs")
rho_cbsa <- cor(state_rr$min_cbsa_dist,state_rr$RR,method= "spearman", use = "complete.obs")
rho_nbdist <- cor(state_rr$nb_dist,state_rr$RR,method = "spearman",use = "complete.obs") #Ignore NAs


fn_state_euclid_dist_plot <- paste0("figs/",scenario,"/state_euclid_dist_plot.jpg") 
state_euclid_dist_plot <- state_rr %>%
  filter(x != y) %>%
  ggplot() +
  geom_point(aes(x=euclid_dist,y=RR),alpha=0.1,size=1,fill='black') +
  geom_smooth(aes(x=euclid_dist,y=RR),color='firebrick',fill='firebrick',method='loess') +
  scale_x_continuous(name="State Centroid Distance (km)",
                     breaks=c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),
                     limits=c(0,2000)) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
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
  labs(subtitle = paste0("Spearman's rho: ",round(rho_euclid,2)))

ggsave(fn_state_euclid_dist_plot,
       plot=state_euclid_dist_plot,
       device = "jpeg",
       dpi = 192,
       width = 6,
       height = 4
)

fn_state_cbsa_dist_plot <- paste0("figs/",scenario,"/state_cbsa_dist_plot.jpg") 
state_cbsa_dist_plot <- state_rr %>%
  filter(x != y) %>%
  ggplot() +
  geom_point(aes(x=min_cbsa_dist,y=RR),alpha=0.1,size=1,fill='black') +
  geom_smooth(aes(x=min_cbsa_dist,y=RR),color='firebrick',fill='firebrick',method='loess') +
  scale_x_continuous(name="State Nearest CBSA Distance (km)",
                     breaks=c(0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000,9500,10000),
                     limits=c(0,2000)) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
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
  labs(subtitle = paste0("Spearman's rho: ",round(rho_cbsa,2)))

ggsave(fn_state_cbsa_dist_plot,
       plot=state_cbsa_dist_plot,
       device = "jpeg",
       dpi = 192,
       width = 6,
       height = 4
)

fn_euclid_cbsa_dist_plot <- paste0("figs/",scenario,"/state_euclid_vs_cbsa.jpg")
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
                     name=expression(RR["identical sequences"]),
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.5),10^(1.5))) +
  geom_hline(yintercept = 1) + #Reference point
  theme_bw()

fn_state_euclid_logdist_plot <- paste0("figs/",scenario,"/state_euclid_logdist_plot.jpg") 
state_euclid_logdist_plot <- state_rr %>%
  ggplot() +
  geom_point(aes(x=log10(euclid_dist+1),y=RR),alpha=0.1,size=1,fill='black') +
  geom_smooth(aes(x=log10(euclid_dist+1),y=RR),color='firebrick',fill='firebrick',method='loess') +
  scale_x_continuous(name=expression(log["10"]("Centroid Distance")),
                     breaks = c(0,1,2,3,4,5),
                     labels = c(0,expression(10^{1}),expression(10^{2}),expression(10^{3}),expression(10^{4}),expression(10^{5})),
                     limits=c(0,4)) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
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
  labs(subtitle = paste0("Spearman's rho: ",round(rho_euclid,2)))

ggsave(fn_state_euclid_logdist_plot,
       plot=state_euclid_logdist_plot,
       device = "jpeg",
       dpi = 192,
       width = 6,
       height = 6
)

fn_state_nb_dist_plot <- paste0("figs/",scenario,"/state_nb_dist_plot.jpg")
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
                     name=expression(RR["identical sequences"]),
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
  labs(subtitle = paste0("Spearman's rho: ",round(rho_nbdist,2)))

ggsave(fn_state_nb_dist_plot,
       plot=state_nb_dist_plot,
       device = "jpeg",
       dpi = 192,
       width = 6,
       height = 4,
       create.dir = TRUE
)

# Create 2x2 patchwork of key distance/travel metrics
# Add same_state indicator
state_rr <- state_rr %>%
  mutate(same_state = ifelse(x == y, "Same state", "Different state"))

# Calculate correlations
rho_cbsa <- cor(state_rr$min_cbsa_dist, state_rr$RR, method = "spearman", use = "complete.obs")
rho_trips <- cor(state_rr$RR_trips, state_rr$RR, method = "spearman", use = "complete.obs")
rho_move <- cor(state_rr$RR_move, state_rr$RR, method = "spearman", use = "complete.obs")

# 1. Neighbor rank plot with ANOVA (miniaturized)
# Perform ANOVA
nb_data <- state_rr %>% filter(x != y, !is.na(nb_dist), nb_dist <= 6)
anova_result <- aov(log10(RR) ~ as.factor(nb_dist), data = nb_data)
anova_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
anova_text <- ifelse(anova_p < 0.001, "p < 0.001",
                     paste0("p = ", formatC(anova_p, digits = 3, format = "f")))

p1 <- nb_data %>%
  ggplot(aes(x = nb_dist, y = RR, group = nb_dist)) +
  geom_boxplot(fill = NA, linewidth = 0.3, outlier.size = 0.5) +
  scale_x_continuous(name = "Neighbor Order",
                     breaks = 1:6,
                     limits = c(0.5, 6.5)) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR[identical]),
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(subtitle = paste0("ANOVA: ", anova_text)) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9))

# 2. CBSA distance plot (miniaturized) - colored by same/different state
p2 <- state_rr %>%
  ggplot(aes(x = min_cbsa_dist, y = RR)) +
  geom_point(aes(color = same_state), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = c("Same state" = "cornflowerblue", "Different state" = "lightcoral"),
                     name = "") +
  scale_x_continuous(name = "CBSA Distance (km)",
                     limits = c(0, 3000)) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR[identical]),
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(subtitle = paste0("ρ = ", round(rho_cbsa, 2))) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# 3. NHTS trips plot (miniaturized) - colored by same/different state
# Get full data extent
trips_range <- range(state_rr$RR_trips[!is.na(state_rr$RR_trips) & state_rr$RR_trips > 0])
p3 <- state_rr %>%
  filter(!is.na(RR_trips)) %>%
  ggplot(aes(x = RR_trips, y = RR)) +
  geom_point(aes(color = same_state), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = c("Same state" = "cornflowerblue", "Different state" = "lightcoral"),
                     name = "") +
  scale_x_continuous(name = "Travel RR (NHTS)",
                     transform = 'log10',
                     breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                     labels = c(expression(10^{-3}), expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}), expression(10^{2}))) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR[identical]),
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(subtitle = paste0("ρ = ", round(rho_trips, 2))) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# 4. SafeGraph mobility plot (miniaturized) - colored by same/different state
# Get full data extent
move_range <- range(state_rr$RR_move[!is.na(state_rr$RR_move) & state_rr$RR_move > 0])
p4 <- state_rr %>%
  filter(!is.na(RR_move)) %>%
  ggplot(aes(x = RR_move, y = RR)) +
  geom_point(aes(color = same_state), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = c("Same state" = "cornflowerblue", "Different state" = "lightcoral"),
                     name = "") +
  scale_x_continuous(name = "Mobility RR (SafeGraph)",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10, 100),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}), expression(10^{2}))) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR[identical]),
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(subtitle = paste0("ρ = ", round(rho_move, 2))) +
  theme(plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# Combine into 2x2 grid with shared legend
combined_dist_plot <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fn_combined_dist <- paste0("figs/", scenario, "/state_distance_poster.jpg")
ggsave(fn_combined_dist,
       plot = combined_dist_plot,
       device = "jpeg",
       dpi = 192,
       width = 10,
       height = 8,
       create.dir = TRUE
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
epsilon <- 1E-4
rho_short_air <- state_rr_air %>%
  filter(flight_length == "Short") %>%
  summarize(rho = cor(RR_air + epsilon, RR + epsilon, method = "spearman", use = "complete.obs")) %>%
  pull(rho)

rho_medium_air <- state_rr_air %>%
  filter(flight_length == "Medium") %>%
  summarize(rho = cor(RR_air + epsilon, RR + epsilon, method = "spearman", use = "complete.obs")) %>%
  pull(rho)

rho_long_air <- state_rr_air %>%
  filter(flight_length == "Long") %>%
  summarize(rho = cor(RR_air + epsilon, RR + epsilon, method = "spearman", use = "complete.obs")) %>%
  pull(rho)

# Create plot for short flights
p_air_short <- state_rr_air %>%
  filter(flight_length == "Short") %>%
  ggplot(aes(x = RR_air + epsilon, y = RR + epsilon)) +
  geom_point(aes(color = same_state), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = c("Same state" = "cornflowerblue", "Different state" = "lightcoral"),
                     name = "") +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR[identical]),
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Short distance (<300 km)",
       subtitle = paste0("ρ = ", round(rho_short_air, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "none")

# Create plot for medium flights
p_air_medium <- state_rr_air %>%
  filter(flight_length == "Medium") %>%
  ggplot(aes(x = RR_air + epsilon, y = RR + epsilon)) +
  geom_point(aes(color = same_state), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = c("Same state" = "cornflowerblue", "Different state" = "lightcoral"),
                     name = "") +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR[identical]),
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Medium distance (300-1800 km)",
       subtitle = paste0("ρ = ", round(rho_medium_air, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "none")

# Create plot for long flights
p_air_long <- state_rr_air %>%
  filter(flight_length == "Long") %>%
  ggplot(aes(x = RR_air + epsilon, y = RR + epsilon)) +
  geom_point(aes(color = same_state), alpha = 0.5, size = 1) +
  geom_smooth(method = 'loess', linewidth = 1.5, se = FALSE, color = alpha("black", 0.7)) +
  scale_color_manual(values = c("Same state" = "cornflowerblue", "Different state" = "lightcoral"),
                     name = "") +
  scale_x_continuous(name = "Air Travel RR",
                     transform = 'log10',
                     breaks = c(0.01, 0.1, 1, 10),
                     labels = c(expression(10^{-2}), expression(10^{-1}),
                               expression(10^{0}), expression(10^{1}))) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR[identical]),
                     breaks = c(1E-1, 1, 1E1),
                     labels = c(expression(10^{-1}), expression(10^{0}), expression(10^{1})),
                     limits = c(10^(-1.5), 10^(1.5))) +
  geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
  theme_classic(base_size = 10) +
  labs(title = "Long distance (>1800 km)",
       subtitle = paste0("ρ = ", round(rho_long_air, 2))) +
  theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "none")

# Combine into 3x1 vertical grid
combined_air_plot <- p_air_short / p_air_medium / p_air_long

fn_combined_air <- paste0("figs/", scenario, "/state_air_distance_poster.jpg")
ggsave(fn_combined_air,
       plot = combined_air_plot,
       device = "jpeg",
       dpi = 192,
       width = 5,
       height = 8,
       create.dir = TRUE
)

message(paste("Air travel distance poster saved to:", fn_combined_air))
