#File: census_div_heatmap.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Makes a heatmap from census divisions RR matrix
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
#Manual ordering for labels
#DIV_ORDER<- c("Pacific","Mountain","West North Central","East North Central","Middle Atlantic",
#"New England","South Atlantic","East South Central","West South Central","Alaska","Hawaii","Pacific Territories","Caribbean Territories")
DIV_ORDER <- c("Far West","Rocky Mountain","Southwest","Southeast","Plains","Great Lakes","Mideast","New England","Mexico","Western Canada","Eastern Canada")

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_census_div.tsv",sep="")
div_rr <- fread(fn_rr)

UB<-2 #Set as bounds for RR and transform to log for display purposes
LB <-1/UB
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}

div_rr$x <- factor(div_rr$x,levels=DIV_ORDER)
div_rr$y <- factor(div_rr$y,levels=DIV_ORDER)

div_heatmap <- div_rr %>% rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
  ggplot(aes(x=x,y=y,fill=fill_RR)) +
  geom_tile() +
  RR_log_grad(LB,UB) +
  theme_minimal() + theme(plot.title=element_text(hjust=0.5)) + 
  geom_text_repel(aes(label=round(RR,digits=2)),color = "black",
                  size=2,
                  bg.color="white",bg.r=0.1,
                  force = 0)+
  labs(x="Regions",y="Regions",title=scenario) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
fn_div_plot <- paste0("figs/",scenario,"/census_divisions_heatmap.jpg") 

ggsave(fn_div_plot,
       plot=div_heatmap,
       device = "jpeg",
       dpi = 192,
       width = 7,
       height = 6
)
