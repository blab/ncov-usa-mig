#File: state_dist_plot.R
#Author(s): Amin Bemanian
#Date: 8/20/24
#Description: Makes a heatmap from state RR matrix
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggsignif)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
state_rr <- fread(fn_rr)  %>%
  mutate(simple_adj = case_when(
    nb_dist == 0 ~ "Within state",
    nb_dist == 1 ~ "Adjacent state",
    TRUE ~ "Non-adjacent state"
  ))  

state_rr$simple_adj<-factor(state_rr$simple_adj,levels=c("Non-adjacent state","Adjacent state","Within state"))

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
              map_signif_level = TRUE,
              color = 'firebrick') +
  geom_boxplot(data = df_adj_boxplot,
               aes(ymin = y05_plot, lower = y25_plot, middle = y50_plot,
                   upper = y75_plot, ymax = y95_plot),
               stat = "identity", fill = NA, width = 0.7) +
  scale_x_discrete(name= '') +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(1E-2,1E-1, 1, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-1.001),10^(1.5))) +
  coord_flip() +
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  ggtitle(scenario)  

ggsave(fn_state_adj_plot,
       plot=state_adj_plot,
       device = "jpeg",
       dpi = 600,
       width = 6,
       height = 6
)

#Note since spearman is rank-order correlation, log transformation does not matter
rho_euclid <- cor(state_rr$euclid_dist,state_rr$RR,method = "spearman")
rho_nbdist <- cor(state_rr$nb_dist,state_rr$RR,method = "spearman",use = "complete.obs") #Ignore NAs


fn_state_euclid_dist_plot <- paste0("figs/",scenario,"/state_euclid_dist_plot.jpg") 
state_euclid_dist_plot <- state_rr %>%
  ggplot() +
  geom_point(aes(x=euclid_dist,y=RR),alpha=0.1,size=1,fill='black') +
  geom_linerange(aes(x=euclid_dist,y=RR,ymin = ci_lb, ymax = ci_ub),alpha=0.1,color='black') +
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
       dpi = 600,
       width = 6,
       height = 6
)

fn_state_euclid_logdist_plot <- paste0("figs/",scenario,"/state_euclid_logdist_plot.jpg") 
state_euclid_logdist_plot <- state_rr %>%
  ggplot() +
  geom_point(aes(x=log10(euclid_dist+1),y=RR),alpha=0.1,size=1,fill='black') +
  geom_linerange(aes(x=euclid_dist,y=RR,ymin = ci_lb, ymax = ci_ub),alpha=0.1,color='black') +
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
       dpi = 600,
       width = 6,
       height = 6
)

fn_state_nb_dist_plot <- paste0("figs/",scenario,"/state_nb_dist_plot.jpg")
nb_jitter <- position_jitter(width=0.2,height=0)
state_nb_dist_plot <- state_rr %>%
  ggplot() +
  #geom_linerange(aes(x=nb_dist,y=RR,ymin = ci_lb, ymax = ci_ub),
   #               position = nb_jitter,
    #              alpha = 0.1,color='black') +
  geom_point(aes(x=nb_dist,y=RR),
             position = nb_jitter,
             alpha = 0.1,fill='black',size=1) +
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
       dpi = 600,
       width = 6,
       height = 6
)
