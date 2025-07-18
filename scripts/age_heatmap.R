#File: age_heatmap.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Makes a heatmap from age RR matrix
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
scenario <- "CAM_1000"

if(scenario ==  "CAM_1000"){
  SCALE_FACTOR <- 2.5
  RR_SIZE <- 0
  AXIS_SIZE <- 8
}else{
  SCALE_FACTOR <- 1
  RR_SIZE <- 1.5
  AXIS_SIZE <- 12
}

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_age_class.tsv",sep="")
age_rr <- fread(fn_rr)

#Force a fill value to be within a certain range for display purposes 
UB<-1.3 #Set as bounds for RR and transform to log for display purposes
LB <- 2-UB #Symmetrical bounds around 1
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}


age_heatmap <- age_rr %>% rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
  ggplot(aes(x=x,y=y,fill=fill_RR)) +
  geom_tile() +
  RR_log_grad(LB,UB) +
  theme_minimal() + theme(plot.title=element_text(hjust=0.5)) + 
  geom_text_repel(aes(label=round(RR,digits=2)),color = "black",
                  size=RR_SIZE,
                  bg.color="white",bg.r=0.1,
                  force = 0)+
  labs(x="Age Groups",y="Age Groups",title=scenario) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = AXIS_SIZE),
        axis.text.y = element_text(size = AXIS_SIZE))

fn_age_plot <- paste0("figs/",scenario,"/age_heatmap",".jpg") 
  
ggsave(fn_age_plot,
       plot=age_heatmap,
       device = "jpeg",
       dpi = 192,
       width = 25,
       height = 20
)