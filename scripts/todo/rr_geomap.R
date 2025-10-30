#File: rr_geomap.R
#Author(s): Amin Bemanian
#Date: 10/2/24
#Description: Geographic map of RRs, will iteratively go over each geo-subunit and make a map for each
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--scale: "state" vs "county" (Later addition if we can do geocoding)
#TODO: Adapt for CAM scenarios, or replace with D3/Shiny and discontinue

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(usmap)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--scale', type = 'character', default = 'state', help = 'Level of geography (state vs county)')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
scale <- args$scale #For now county will not be implemented
scenario <- "USA"

dir.create(paste("figs/",scenario,"/maps",sep=""))

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
geo_rr <- fread(fn_rr)

UB<-1.2 #Set as bounds for RR and transform to log for display purposes
LB <- 2-UB #Symmetrical bounds around 1
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}

geo_list <- unique(geo_rr$x)
for(g in geo_list){
  g_mat <- filter(geo_rr,x==g) %>%
    select(y,RR) %>%
    rename(state = y) %>%
    rowwise %>% mutate(values = fill_bound(RR))
  fn_map <- paste0("figs/",scenario,"/maps/rr_map_",g,".jpg")
  rr_map <- plot_usmap(regions = scale,data=g_mat) +
    scale_fill_gradient2(name="RR",
                         high = "#D67C34",
                         low = "#4C90C0",
                         limits=c(log10(LB),log10(UB)),
                         breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                         labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB))+
    theme(plot.title=element_text(hjust=0.5)) +
    labs(title = g)
  ggsave(fn_map,
         plot=rr_map,
         device = "jpeg",
         dpi = 600,
         width = 6,
         height = 4
  )
}
