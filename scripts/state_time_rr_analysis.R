library(tidyverse)
library(data.table)
library(dbplyr)
library(duckdb)
library(ggplot2)
library(tictoc)
library(scales)

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")

select <- dplyr::select
STATE_DISTANCES <- fread("data/nb_dist_states.tsv")
SCENARIO <- "USA_mini"

STATE_ORDER<- c("Alaska","Hawaii","Washington","Oregon","California",
                "Nevada","Idaho","Montana","Wyoming","Utah","Colorado","Arizona","New Mexico",
                "North Dakota","South Dakota","Nebraska","Kansas","Missouri","Iowa","Minnesota",
                "Wisconsin","Illinois","Indiana","Michigan","Ohio",
                "Pennsylvania","New York","New Jersey",
                "Delaware","Maryland","District of Columbia","Virginia",
                "Connecticut","Rhode Island","Massachusetts","Vermont","New Hampshire","Maine",
                "Tennessee","Kentucky","West Virginia",
                "North Carolina","South Carolina","Georgia","Florida","Alabama","Mississippi",
                "Arkansas","Louisiana","Texas","Oklahoma",
                "Guam","American Samoa","Puerto Rico","Virgin Islands")


fn_db <- paste0("db_files/db_",SCENARIO,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)
df_meta <- tbl(con,"metadata")
df_pairs <- tbl(con,"pairs_time")

strain_counts <- df_pairs %>% 
  group_by(strain_1) %>% 
  summarize(
    n = n(),
    Nextstrain_clade = first(Nextstrain_clade),
  ) %>% 
  arrange(desc(n)) %>%
  collect()

hist_max <- 2.5E4
strain_hist <- ggplot(data = strain_counts,aes(x = n)) +
  geom_histogram(fill="lightblue",color="black",bins=20) +
  labs(x = "Sequence frequency in pair set (log-scaled)",
       y = "Number of sequences") +
  stat_ecdf(
    aes(y = ..y.. * hist_max),
    geom = "line",
    color = "red",
    size = 1,
    lty =2
  ) +
  geom_line(y = 0.90 * hist_max, color = "black", lty = 2) +
  scale_x_log10(labels=label_comma()) +
  scale_y_continuous(
    labels=label_comma(),
    sec.axis = sec_axis(
      trans = ~ . * (100 / hist_max),  # Rescale back to 0–100%
      name = "Cumulative %",
      breaks=seq(0,100,by=20)
    )
  ) +
  labs(title="All clades") +
  theme_bw() 


ggplot(data = filter(strain_counts,Nextstrain_clade == "21K (Omicron)"),aes(x = n)) +
  geom_histogram(fill="lightblue",color="black",bins=20) +
  labs(x = "Sequence frequency in pair set (log-scaled)",
       y = "Number of sequences") +
  stat_ecdf(
    aes(y = ..y.. * hist_max),
    geom = "line",
    color = "red",
    size = 1,
    lty =2
  ) +
  geom_line(y = 0.90 * hist_max, color = "black", lty = 2) +
  scale_x_log10(labels=label_comma()) +
  scale_y_continuous(
    labels=label_comma(),
    sec.axis = sec_axis(
      trans = ~ . * (100 / hist_max),  # Rescale back to 0–100%
      name = "Cumulative %",
      breaks=seq(0,100,by=20)
    )
  ) +
  labs(title="Clade 21K (Early Omicron)") +
  theme_bw() 


hCut <- 500

keep_strains <- strain_counts %>% 
  filter(n < hCut) %>% 
  select(strain_1) %>%
  rename(strains = strain_1)

#Duplicate DB for testing purposes
#WILL DELETE AFTER THIS TESTING ON APR-17 IS DONE
#I AM NOT PROUD BUT THIS IS WHAT IT IS
df_meta <- tbl(con,"metadata")
df_pairs <- tbl(con,"pairs_time")
filtered_pairs <- df_pairs %>% collect %>%
  filter(strain_1 %in% keep_strains$strains) %>%
  filter(strain_2 %in% keep_strains$strains)
dbWriteTable(con,"pairs",filtered_pairs,overwrite=TRUE)
dbWriteTable(con,"pairs_time",filtered_pairs,overwrite=TRUE)

state_rr_all <- con %>%
  bind_pairs_exp("division") %>%
  calculate_rr_matrix()

state_rr_all$x <- factor(state_rr_all$x,levels=STATE_ORDER)
state_rr_all$y <- factor(state_rr_all$y,levels=STATE_ORDER)

#Generalize upper and lowerbound definition
TIME_LB_YEARS <- seq(as.Date("2020-01-01"),as.Date("2025-01-01"),by="6 mo")
TIME_UB_YEARS <- TIME_LB_YEARS[2:length(TIME_LB_YEARS)]-1
TIME_LB_YEARS <- TIME_LB_YEARS[1:length(TIME_LB_YEARS)-1] 

state_rr_year <- NULL

tic("Snapshot analysis")
for(i in 1:length(TIME_LB_YEARS)){
  if(i == 1){
    state_rr_year <- con %>% 
      bind_pairs_exp("division",time_bounds = c(TIME_LB_YEARS[i],TIME_UB_YEARS[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(year_date = TIME_LB_YEARS[i])
  } else{
    rr_year <- con %>% 
      bind_pairs_exp("division",time_bounds = c(TIME_LB_YEARS[i],TIME_UB_YEARS[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(year_date = TIME_LB_YEARS[i])
    state_rr_year <- bind_rows(state_rr_year,rr_year)
  } 
}

#Make sure there are no missing entries and if there are add NAs
all_combinations <- expand.grid(
  x = unique(state_rr_year$x),
  y = unique(state_rr_year$y),
  year_date = unique(state_rr_year$year_date)
)
state_rr_year <- full_join(all_combinations, state_rr_year, by = c("x", "y", "year_date"))
toc()


##Heatmaps
UB<-2 #Set as bounds for RR and transform to log for display purposes
LB <- 0.5 #Symmetrical bounds around 1
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}

SCALE_FACTOR <- 1
RR_SIZE <- 1.5
AXIS_SIZE <- 12

hmap_state_rr_all <- state_rr_all %>%
  rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
  ggplot(aes(x=x,y=y,fill=fill_RR)) +
  geom_tile() +
  scale_fill_gradient2(name="RR",
                       high = "#D67C34",
                       low = "#4C90C0",
                       na.value = "black",
                       limits=c(log10(LB),log10(UB)),
                       breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                       labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB)) +
  theme_minimal() + 
  theme(plot.title=element_text(hjust=0.5,size = 16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = AXIS_SIZE),
        axis.text.y = element_text(size = AXIS_SIZE)) +
  labs(title=paste0("All years - Hyper Threshold: ",hCut))


hmap_state_rr_time <- state_rr_year %>%
  mutate(x = factor(x,levels=STATE_ORDER)) %>%
  mutate(y = factor(y,levels=STATE_ORDER)) %>%
  rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
  ggplot(aes(x=x,y=y,fill=fill_RR)) +
  facet_wrap(facets=vars(year_date),ncol=5,dir='v') +
  geom_tile() +
  scale_fill_gradient2(name="RR",
                       high = "#D67C34",
                       low = "#4C90C0",
                       na.value = "black",
                       limits=c(log10(LB),log10(UB)),
                       breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                       labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB)) +
  theme_minimal() + 
  theme(plot.title=element_text(hjust=0.5,size = 16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16)) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = AXIS_SIZE),
        axis.text.y = element_text(size = AXIS_SIZE)) +
  labs(title=paste0("Hyper threshold: ",hCut))

hmap_state_rr_time
break

# Make date sequences for the time series analysis
lb_start_date <- as.Date("2020-01-01")
lb_end_date <- as.Date("2024-10-08")
lb_vector <- seq(from = lb_start_date, to = lb_end_date, by = "4 weeks")

tic("Series anaylsis")
ub_start_date <- as.Date("2020-03-25")
ub_end_date <- as.Date("2024-12-31")
ub_vector <- seq(from = ub_start_date, to = ub_end_date, by = "4 weeks")

state_rr_series <- NULL
for(i in 1:length(lb_vector)){
  if(i == 1){
    state_rr_series <- con %>% 
      bind_pairs_exp("division",time_bounds = c(lb_vector[i],ub_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(ub_vector[i],lb_vector[i])))
  } else{
    rr_series <- con %>% 
      bind_pairs_exp("division",time_bounds = c(lb_vector[i],ub_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(ub_vector[i],lb_vector[i])))
    state_rr_series <- bind_rows(state_rr_series,rr_series)
  } 
}
toc()

ggplot(data=(state_rr_series |> filter(x==y)),
       aes(x=date,y=RR,color=x)) +
  geom_line() +
  geom_point(aes(fill=x)) +
  geom_line(y=0,lty=2,size=1,color="black") +
  facet_wrap(vars(x),ncol=3) + 
  scale_x_date(date_breaks = "6 months") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size= 11)) +
  theme(axis.text.y = element_text(size=11)) +
  labs(x = "Calendar Date", y = "RR",title="Monthly Same State Enrichment") +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(legend.position = "none")

#Function that makes a 2x2 set of plots looking at self-enrichment vs pair-enrichment
#s1 and s2 are state names to be compared
state_pair_plot <- function(s1,s2){
  plot <- ggplot(data=(state_rr_series |> 
                         filter(x %in% c(s1,s2) & y %in% c(s1,s2))),
                 aes(x=date,y=RR,color=x)) +
    geom_line() +
    geom_point(aes(fill=x)) +
    geom_line(y=0,lty=2,size=1,color="black") +
    facet_grid(rows=vars(x),cols=vars(y)) + 
    scale_x_date(date_breaks = "6 months") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size= 11)) +
    theme(axis.text.y = element_text(size=11)) +
    labs(x = "Calendar Date", y = "RR",
         title=paste0("Monthly Enrichment - ",s1," & ",s2)) +
    theme(axis.title.x = element_text(size=14)) +
    theme(axis.title.y = element_text(size=14)) +
    theme(legend.position = "none")
  return(plot)
}

state_pair_plot("Montana", "Wyoming")

# Early Normalization Code - Will Put on Hiatus While Cécile works on it
# all_combinations <- expand.grid(
#   x = unique(state_rr_series$x),
#   y = unique(state_rr_series$y),
#   date = unique(state_rr_series$date)
# )
# 
# state_rr_series <- full_join(all_combinations, state_rr_series, by = c("x", "y", "date"))
# 
# RR_WA_WA <- state_rr_series |> filter(x=="Washington",y=="Washington") 
# RR_OR_OR <- state_rr_series |> filter(x=="Oregon",y=="Oregon") 
# RR_WA_OR <- state_rr_series |> filter(x=="Washington",y=="Oregon")
# 
# nRR_WA_OR <- data_frame(
#   date = RR_WA_OR$date,
#   RR = RR_WA_OR$RR/sqrt(RR_WA_WA$RR*RR_OR_OR$RR) 
# )
# 
# ggplot(data=(nRR_WA_OR),
#        aes(x=date,y=RR)) +
#   geom_line() +
#   geom_point() +
#   scale_x_date(date_breaks = "6 months") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size= 11)) +
#   theme(axis.text.y = element_text(size=11)) +
#   labs(x = "Calendar Date", y = "Normalzied RR",title="Monthly Enrichment - Washington & Oregon") +
#   theme(axis.title.x = element_text(size=14)) +
#   theme(axis.title.y = element_text(size=14))
# 
# RR_ID_ID <- state_rr_series |> filter(x=="Idaho",y=="Idaho") 
# RR_MT_MT <- state_rr_series |> filter(x=="Montana",y=="Montana") 
# RR_ID_MT <- state_rr_series |> filter(x=="Idaho",y=="Montana")
# 
# nRR_ID_MT <- data_frame(
#   date = RR_ID_MT$date,
#   RR = RR_ID_MT$RR/sqrt(RR_ID_ID$RR*RR_MT_MT$RR) 
# )
# 
# 
# ggplot(data=(nRR_ID_MT),
#        aes(x=date,y=RR)) +
#   geom_line() +
#   geom_point() +
#   scale_x_date(date_breaks = "6 months") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size= 11)) +
#   theme(axis.text.y = element_text(size=11)) +
#   labs(x = "Calendar Date", y = "Normalzied RR",title="Monthly Enrichment - Idaho & Montana") +
#   theme(axis.title.x = element_text(size=14)) +
#   theme(axis.title.y = element_text(size=14))

fn_out_path <- paste0("results/",SCENARIO,"/time_state/")
dir.create(file.path(fn_out_path),showWarnings = FALSE)
write_tsv(state_rr_all,file=paste0(fn_out_path,"df_state_rr_all.tsv"))
write_tsv(state_rr_series,file=paste0(fn_out_path,"df_state_rr_series.tsv"))
write_tsv(state_rr_year,file=paste0(fn_out_path,"df_state_rr_snap.tsv"))

dbDisconnect(con)
