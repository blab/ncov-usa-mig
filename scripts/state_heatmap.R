#File: state_heatmap.R
#Author(s): Amin Bemanian
#Date: 8/15/24
#Description: Makes a heatmap from state RR matrix
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
STATE_ORDER<- c("Alaska","Hawaii","Washington","Oregon","California",
"Nevada","Idaho","Montana","Wyoming","Utah","Colorado","Arizona","New Mexico",
"North Dakota","South Dakota","Nebraska","Kansas","Missouri","Iowa","Minnesota",
"Wisconsin","Illinois","Indiana","Michigan","Ohio",
"Pennsylvania","New York","New Jersey",
"Delaware","Maryland","District of Columbia","Virginia",
"Connecticut","Rhode Island","Massachusetts","Vermont","New Hampshire","Maine",
"Tennessee","Kentucky","West Virginia",
"North Carolina","South Carolina","Georgia","Florida","Alabama","Mississippi",
"Arkansas","Louisiana","Texas","Oklahoma","Mexico",
"British Columbia","Alberta","Saskatchewan","Manitoba",
"Ontario","Quebec","Newfoundland and Labrador","New Brunswick","Prince Edward Island", "Nova Scotia")

if(scenario ==  "USA"){
  SCALE_FACTOR <- 2.5
  RR_SIZE <- 0
  AXIS_SIZE <- 8
}else{
  SCALE_FACTOR <- 2.5
  RR_SIZE <- 0
  AXIS_SIZE <- 8
}

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
state_rr <- fread(fn_rr)
print(unique(state_rr$x))


UB<-1.2 #Set as bounds for RR and transform to log for display purposes
LB <- 2-UB #Symmetrical bounds around 1
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}

print(paste("The range of RRs is:",range(state_rr$RR)))

state_rr$x <- factor(state_rr$x,levels=STATE_ORDER)
state_rr$y <- factor(state_rr$y,levels=STATE_ORDER)

#TEMPORARY FOR PIDS GRANT PROPOSAL DELETE IN FINAL VERSION
#excl_states <- c("Minnesota","Wisconsin","Illinois","Indiana","Michigan","Ohio","Kentucky","West Virginia")
#state_rr <- state_rr |> filter(!(x %in% excl_states)) |> filter(!(y %in% excl_states))

state_heatmap <- state_rr %>% rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
  ggplot(aes(x=x,y=y,fill=fill_RR)) +
  geom_tile() +
  scale_fill_gradient2(name="RR",
                       high = "#D67C34",
                       low = "#4C90C0",
                       limits=c(log10(LB),log10(UB)),
                       breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                       labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB)) +
  theme_minimal() + 
  theme(plot.title=element_text(hjust=0.5,size = 24,face="bold"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + 
  
  #geom_text_repel(aes(label=round(RR,digits=2)),color = "black", #Re-enable if you want RR to be seen
  #                size=RR_SIZE,
  #                bg.color="white",bg.r=0.1,
  #                force = 0)+
  labs(x="Divisions",y="Divisions",title=scenario) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = AXIS_SIZE),
        axis.text.y = element_text(size = AXIS_SIZE))
fn_state_plot <- paste0("figs/",scenario,"/state_heatmap",".jpg") 

ggsave(fn_state_plot,
       plot=state_heatmap,
       device = "jpeg",
       dpi = 200,
       width = 7 * SCALE_FACTOR,
       height = 6 * SCALE_FACTOR,
       create.dir = TRUE
)
