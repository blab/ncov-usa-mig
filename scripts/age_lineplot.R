#File: age_lineplot.R
#Author(s): Amin Bemanian
#Date: 10/7/24
#Description: Makes a lineplot of the relative risk for each age group vs other age groups
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_age_class.tsv",sep="")
age_rr <- fread(fn_rr)
fn_rr_adj <- paste("results/",scenario,"/df_RR_by_age_adj.tsv",sep="")
age_rr_adj <- fread(fn_rr_adj)

my_pal_age <- viridis_pal(option = 'D', direction = -1)(
  length(unique(age_rr$x)) #Pass the length
)


#x is 1st age group, y is 2nd age group
age_plot_all <- age_rr %>% 
  ggplot(aes(x = x, colour = as.factor(y))) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = y)) +
  #geom_linerange(aes(ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 4, scales = 'free_y') +
  scale_x_discrete(name = 'Age Group') +
  scale_y_continuous(name = expression(RR["identical sequences"])) +
  scale_colour_manual(values = my_pal_age) + 
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario," - All Pairs")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5))

fn_age_plot_all <- paste0("figs/",scenario,"/age_lineplot_all.jpg") 
ggsave(fn_age_plot_all,
       plot=age_plot_all,
       device = "jpeg",
       dpi = 600,
       width = 24,
       height = 24
)

#Do RRs stratified by adjacency plot (single plot)
age_plot_by_adj <- age_rr_adj %>% 
  ggplot(aes(x = x, colour = as.factor(adj_status))) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = adj_status)) +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 2, scales = 'free_y') +
  scale_x_discrete(name = 'Age Group') +
  scale_y_continuous(name = expression(RR["identical sequences"])) +
  scale_color_brewer(type="qual", name = "State Adjacency",palette = "Paired") + 
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario," - By Adjacency")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        plot.title=element_text(hjust=0.5)) 
  
fn_age_plot_by_adj <- paste0("figs/",scenario,"/age_lineplot_by_adj.jpg") 
ggsave(fn_age_plot_by_adj,
       plot=age_plot_by_adj,
       device = "jpeg",
       dpi = 600,
       width = 8,
       height = 10
)

#Now breakdown adjacency stratified RRs into 3 plots
age_plot_within <- age_rr_adj %>% filter(adj_status == "Within") %>% 
  ggplot(aes(x = x, colour = as.factor(y))) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = y)) +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 2, scales = 'free_y') +
  scale_x_discrete(name = 'Age Group') +
  scale_y_continuous(name = expression(RR["identical sequences"])) +
  scale_colour_manual(values = my_pal_age) + 
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario," - Within Same State")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5))

fn_age_plot_within <- paste0("figs/",scenario,"/age_lineplot_within.jpg") 
ggsave(fn_age_plot_within,
       plot=age_plot_within,
       device = "jpeg",
       dpi = 600,
       width = 7,
       height = 10
)

age_plot_adjacent <- age_rr_adj %>% filter(adj_status == "Adjacent") %>% 
  ggplot(aes(x = x, colour = as.factor(y))) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = y)) +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 2, scales = 'free_y') +
  scale_x_discrete(name = 'Age Group') +
  scale_y_continuous(name = expression(RR["identical sequences"])) +
  scale_colour_manual(values = my_pal_age) + 
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario," - Adjacent States")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5))

fn_age_plot_adjacent <- paste0("figs/",scenario,"/age_lineplot_adjacent.jpg") 
ggsave(fn_age_plot_adjacent,
       plot=age_plot_adjacent,
       device = "jpeg",
       dpi = 600,
       width = 7,
       height = 10
)

age_plot_nonadj <- age_rr_adj %>% filter(adj_status == "Non-adjacent") %>% 
  ggplot(aes(x = x, colour = as.factor(y))) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = y)) +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 2, scales = 'free_y') +
  scale_x_discrete(name = 'Age Group') +
  scale_y_continuous(name = expression(RR["identical sequences"])) +
  scale_colour_manual(values = my_pal_age) + 
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario," - Non-adjacent States")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5))

fn_age_plot_nonadj <- paste0("figs/",scenario,"/age_lineplot_nonadj.jpg") 
ggsave(fn_age_plot_nonadj,
       plot=age_plot_nonadj,
       device = "jpeg",
       dpi = 600,
       width = 7,
       height = 10
)
