#File: age_histo.R
#Author(s): Amin Bemanian
#Date: 10/7/24
#Description: Makes a histogram of the age distribution for both the individuals and the pairs
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

my_pal_age <- viridis_pal(option = 'D', direction = -1)(
  length(unique(age_rr$x)) #Pass the length
)

fn_meta <- paste("results/",scenario,"/metadata_age_clean.tsv",sep="")
df_meta <- fread(fn_meta)

my_pal_age <- viridis_pal(option = 'D', direction = -1)(
  length(unique(age_rr$x)) #Pass the length
)

age_pairs_dist <- age_rr %>% 
  group_by(x) %>%
  rename(age_class = x) %>%
  summarise(n = sum(n_pairs)) %>%
  mutate(perc = 100*n/sum(n))

age_seq_dist <- df_meta %>%
  count(age_class) %>%
  mutate(perc = 100*n/sum(n))

plot_seq_hist <- ggplot(data = age_seq_dist,aes(x=age_class,y=n)) +
  geom_col(aes(fill=age_class)) +
  geom_text(aes(label=paste0(round(perc,1),"%"),),hjust=0) +
  scale_fill_manual(values = my_pal_age) + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5)) +
  labs(
    title = "Age Distribution of Individuals",
    x = "Age Group",
    y = "Count",
  ) +
  coord_flip(clip = "off") 

fn_seq_hist <- paste0("figs/",scenario,"/age_seq_hist.jpg")
ggsave(fn_seq_hist,plot_seq_hist,width=9,height = 9, dpi = 600)

plot_pairs_hist <- ggplot(data = age_pairs_dist,aes(x=age_class,y=n)) +
  geom_col(aes(fill=age_class)) +
  geom_text(aes(label=paste0(round(perc,1),"%"),),hjust=0) +
  scale_fill_manual(values = my_pal_age) + 
  theme_classic() +
  theme(legend.position = "none",
        plot.title=element_text(hjust=0.5)) +
  labs(
    title = "Age Distribution of Identical Pairs",
    x = "Age Group",
    y = "Count",
  ) +
  coord_flip(clip = "off")
fn_pairs_hist <- paste0("figs/",scenario,"/age_pairs_hist.jpg")
ggsave(fn_pairs_hist,plot_pairs_hist,width=9,height = 9, dpi = 600)
