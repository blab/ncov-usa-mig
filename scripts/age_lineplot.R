#File: age_lineplot.R
#Author(s): Amin Bemanian
#Date: 07/07/25
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
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000", help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

AGE_UB <- "80y" #Censor older ages group

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_age_class.tsv",sep="")
age_rr <- fread(fn_rr) %>%
  mutate(age_censor = (x > AGE_UB | y > AGE_UB))
fn_age_state_rr <- paste("results/",scenario,"/df_RR_by_age_state.tsv",sep="")
age_state_rr <- fread(fn_age_state_rr) %>%
  mutate(age_censor = (x > AGE_UB | y > AGE_UB))

NTH_BREAK <- 5
age_levels <- age_rr %>% 
  pull(x) %>% 
  factor() %>% 
  levels()

age_breaks <- age_levels[seq(1,length(age_levels),by=NTH_BREAK)] 

my_pal_age <- viridis_pal(option = 'D', direction = -1)(
  length(unique(age_rr$x)) #Pass the length
)

age_group_to_numeric <- function(a){
  gsub("y","",a) |> as.numeric() #Drop the y and make numeric, only works if single age categories
}

#x is 1st age group, y is 2nd age group
age_plot_all <- age_rr %>%
  filter(!age_censor) %>%
  ggplot(aes(x = x, colour = as.factor(y))) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = y)) +
  #geom_linerange(aes(ymin = ci_lb, ymax = ci_ub)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 10, scales = 'free_y') +
  scale_x_discrete(name = 'Age Group', breaks = age_breaks) +
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
       dpi = 192,
       width = 50,
       height = 40,
       limitsize = FALSE
)

# Abridged version: focus on key age groups
TARGET_AGES <- c("05y", "20y", "40y", "70y")

age_plot_abridged <- age_rr %>%
  filter(!age_censor) %>%
  filter(y %in% TARGET_AGES) %>%
  ggplot(aes(x = x,color = y)) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR,group = 1)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 2, scales = 'free_y') +
  scale_x_discrete(name = 'Age Group', breaks = age_breaks) +
  scale_y_continuous(name = expression(RR["identical sequences"]),transform = "log10") +
  theme_classic() +
  ggtitle(paste0("Age-wise RR for ",scenario," - Key Age Groups")) +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        legend.position = 'bottom',
        plot.title=element_text(hjust=0.5))

fn_age_plot_abridged <- paste0("figs/",scenario,"/age_lineplot_abridged.jpg")
ggsave(fn_age_plot_abridged,
       plot=age_plot_abridged,
       device = "jpeg",
       dpi = 300,
       width = 8,
       height = 5
)

age_state_rr <- age_state_rr %>%
  mutate(
    state_region_flag = case_when(
      sameState ~ "Same State",
      !sameState & sameRegion ~ "Same Region, Different States",
      !sameRegion ~ "Different Region",
      .default = "Invalid Region/State"
    ) %>% factor(levels = c("Same State","Same Region, Different States","Different Region"))
  ) 

#Do RRs stratified by adjacency plot (single plot)
age_plot_by_state <- age_state_rr  %>% 
  filter(!age_censor) %>%
  ggplot(aes(x = x, colour = state_region_flag)) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = state_region_flag)) +
  geom_hline(yintercept = 1) +
  facet_wrap(. ~ y, ncol = 10, scales = 'free_y') + #y here refers to age group, not y-axis
  scale_x_discrete(name = 'Age Group', breaks = age_breaks) +
  scale_y_continuous(name = expression(RR["identical sequences"]),transform = "log10") +
  scale_color_brewer(type="qual", name = "Within State",palette = "Dark2") + 
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
       dpi = 192,
       width = 50,
       height = 40,
       limitsize = FALSE
)

##Do a single plot only looking at same age group RR 
age_same_plot_by_state  <- age_state_rr  %>%
  filter(x==y) %>% #Only look at RR within the same 
  filter(!age_censor) %>%
  ggplot(aes(x = x, colour = state_region_flag)) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = state_region_flag)) +
  #geom_linerange(aes(ymin = ci_lb, ymax = ci_ub,group = sameState)) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = 'Age Group', breaks = age_breaks) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(1E-1, 8E-1, 1, 5, 2, 1E1, 1E2),
                     labels = c(0.1,0.8,1,5,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.1),10^(1.8))) +
  scale_color_brewer(type="qual", name = "State/Region Status",palette = "Dark2") + 
  theme_classic() +
  ggtitle(paste0("Same Age RR")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        plot.title=element_text(hjust=0.5),
        legend.position = "bottom"
      )

fn_age_same_state_plot <- paste0("figs/",scenario,"/age_same_lineplot_by_state.jpg")
ggsave(fn_age_same_state_plot,
       plot=age_same_plot_by_state,
       device = "jpeg",
       dpi = 300,
       width = 7,
       height = 4
)

## Boxplot of same age RR by region status
age_same_boxplot_by_state <- age_state_rr %>%
  filter(x==y) %>% #Only look at RR within the same age
  filter(!age_censor) %>%
  ggplot(aes(x = state_region_flag, y = RR, fill = state_region_flag)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(transform = 'log10',
                     name = expression(RR["identical sequences"])) +
  scale_x_discrete(name = "") +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  theme_classic() +
  ggtitle("Same Age RR") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

fn_age_same_boxplot <- paste0("figs/",scenario,"/age_same_boxplot_by_state.jpg")
ggsave(fn_age_same_boxplot,
       plot = age_same_boxplot_by_state,
       device = "jpeg",
       dpi = 300,
       width = 3,
       height = 5
)
