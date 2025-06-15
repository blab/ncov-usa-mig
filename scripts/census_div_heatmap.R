#File: census_div_heatmap.R
#Author(s): Amin Bemanian
#Date: 8/15/24
#Description: Makes a heatmap from census divisions RR matrix
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)


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

UB<-1.2 #Set as bounds for RR and transform to log for display purposes
LB <- 2-UB #Symmetrical bounds around 1
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}

print(paste("The range of RRs is:",range(div_rr$RR)))

div_rr$x <- factor(div_rr$x,levels=DIV_ORDER)
div_rr$y <- factor(div_rr$y,levels=DIV_ORDER)

div_heatmap <- div_rr %>% rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
  ggplot(aes(x=x,y=y,fill=fill_RR)) +
  geom_tile() +
  scale_fill_gradient2(name="RR",
                       high = "#D67C34",
                       low = "#4C90C0",
                       limits=c(log10(LB),log10(UB)),
                       breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                       labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB)) +
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
       dpi = 600,
       width = 7,
       height = 6
)
