#File: age_lineplot.R
#Author(s): Amin Bemanian
#Date: 11/12/24
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
scenario <- "CAM_1000"

AGE_UB <- "79y" #Censor older ages group

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_age_class.tsv",sep="")
age_rr <- fread(fn_rr) %>%
  mutate(age_censor = (x > AGE_UB | y > AGE_UB))
fn_age_state_rr <- paste("results/",scenario,"/df_RR_by_age_state.tsv",sep="")
age_state_rr <- fread(fn_age_state_rr) %>%
  mutate(age_censor = (x > AGE_UB | y > AGE_UB))

my_pal_age <- viridis_pal(option = 'D', direction = -1)(
  length(unique(age_rr$x)) #Pass the length
)

age_group_to_numeric <- function(a){
  gsub("y","",a) |> as.numeric() #Drop the y and make numeric, only works if single age categories
}

#USA is the one scenario with single age categories
if(scenario == "USA"){
  age_rr <- age_rr |>
    mutate(x = age_group_to_numeric(x)) |>
    mutate(y = age_group_to_numeric(y))
}

#x is 1st age group, y is 2nd age group
age_plot_all <- age_rr %>%
  filter(!age_censor) %>%
  ggplot(aes(x = x, colour = as.factor(y))) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = y)) +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 10, scales = 'free_y') +
  scale_x_continuous(name = 'Age Group') +
  scale_y_continuous(name = expression(RR["identical sequences"])) +
  scale_colour_manual(values = my_pal_age) + 
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario," - All Pairs")) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5))

fn_age_plot_all <- paste0("figs/",scenario,"/age_lineplot_all_censor.jpg") 
ggsave(fn_age_plot_all,
       plot=age_plot_all,
       device = "jpeg",
       dpi = 500,
       width = 50,
       height = 40,
       limitsize = FALSE
)

#Do RRs stratified by adjacency plot (single plot)
age_plot_by_state <- age_state_rr  %>% 
  ggplot(aes(x = x, colour = sameState)) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = sameState)) +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub,group = sameState)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 10, scales = 'free_y') + #y here refers to age group, not y-axis
  scale_x_categorical(name = 'Age Group') +
  scale_y_continuous(name = expression(RR["identical sequences"])) +
  scale_color_brewer(type="qual", name = "Within State",palette = "Paired") + 
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario)) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        plot.title=element_text(hjust=0.5)) 
  
fn_age_state_plot <- paste0("figs/",scenario,"/age_lineplot_by_state.jpg") 
ggsave(fn_age_state_plot,
       plot=age_plot_by_state,
       device = "jpeg",
       dpi = 500,
       width = 50,
       height = 40,
       limitsize = FALSE
)

##Do a single plot only looking at same age group RR 
age_same_plot_by_state  <- age_state_rr  %>%
  filter(x==y) %>% #Only look at RR within the same 
  ggplot(aes(x = x, colour = sameState)) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = sameState)) +
  geom_linerange(aes(ymin = ci_lb, ymax = ci_ub,group = sameState)) +
  geom_hline(yintercept = 1) +
  scale_x_continuous(name = 'Age Group') +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(1E-1, 8E-1, 1, 5, 2, 1E1, 1E2),
                     labels = c(0.1,0.8,1,5,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.05),10^(1.2))) +
  scale_color_brewer(type="qual", name = "Within State",palette = "Paired") + 
  theme_classic() +
  ggtitle(paste0("Same Age RR for ",scenario)) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        plot.title=element_text(hjust=0.5))

fn_age_same_state_plot <- paste0("figs/",scenario,"/age_same_lineplot_by_state.jpg") 
ggsave(fn_age_same_state_plot,
       plot=age_same_plot_by_state,
       device = "jpeg",
       dpi = 600,
       width = 8,
       height = 6
)
