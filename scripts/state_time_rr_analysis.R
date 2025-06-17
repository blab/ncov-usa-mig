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
SCENARIO <- "CAM_1000"

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

normalized_state_rr <- function(df_rr){
  # First: filter to the diagonal (RR(x,x)) entries
  rr_diag <- df_rr %>%
    filter(x == y) %>%
    select(date, state = x, RR_diag = RR)
  # Join to get RR(x,x) and RR(y,y) for each row
  df_rr <- df_rr %>%
    left_join(rr_diag, by = c("date", "x" = "state")) %>%
    rename(RR_xx = RR_diag) %>%
    left_join(rr_diag, by = c("date", "y" = "state")) %>%
    rename(RR_yy = RR_diag) %>%
    mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = TRUE))
  return(df_rr)
}

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

state_rr_all <- con %>%
  bind_pairs_exp("division") %>%
  calculate_rr_matrix()

state_rr_all$x <- factor(state_rr_all$x,levels=STATE_ORDER)
state_rr_all$y <- factor(state_rr_all$y,levels=STATE_ORDER)

#Generalize upper and lowerbound definition
TIME_LB_YEARS <- seq(as.Date("2020-01-01"),as.Date("2025-01-01"),by="6 mo")
TIME_UB_YEARS <- TIME_LB_YEARS[2:length(TIME_LB_YEARS)]-1
TIME_LB_YEARS <- TIME_LB_YEARS[1:length(TIME_LB_YEARS)-1] 

state_rr_snap <- NULL

tic("Snapshot analysis")
for(i in 1:length(TIME_LB_YEARS)){
  if(i == 1){
    state_rr_snap <- con %>% 
      bind_pairs_exp("division",time_bounds = c(TIME_LB_YEARS[i],TIME_UB_YEARS[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = TIME_LB_YEARS[i])
  } else{
    rr_snap <- con %>% 
      bind_pairs_exp("division",time_bounds = c(TIME_LB_YEARS[i],TIME_UB_YEARS[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = TIME_LB_YEARS[i])
    state_rr_snap <- bind_rows(state_rr_snap,rr_snap)
  } 
}

#Make sure there are no missing entries and if there are add NAs
all_combinations <- expand.grid(
  x = unique(state_rr_snap$x),
  y = unique(state_rr_snap$y),
  date = unique(state_rr_snap$date)
)
state_rr_snap <- full_join(all_combinations, state_rr_snap, by = c("x", "y", "date"))
toc()

fn_out_path <- paste0("results/",SCENARIO,"/time_state/")
dir.create(file.path(fn_out_path),showWarnings = FALSE)
write_tsv(state_rr_all,file=paste0(fn_out_path,"df_state_rr_all.tsv"))
state_rr_snap <- normalized_state_rr(state_rr_snap)
write_tsv(state_rr_snap,file=paste0(fn_out_path,"df_state_rr_snap.tsv"))

##Heatmaps
UB<-2 #Set as bounds for RR and transform to log for display purposes
LB <- 0.5 #Symmetrical bounds around 1
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}

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

state_rr_series <- normalized_state_rr(state_rr_series)
write_tsv(state_rr_series,file=paste0(fn_out_path,"df_state_rr_series.tsv"))
dbDisconnect(con)

outlier_nRR <- 1.5 #Maximum threshold
ggplot(data=(state_rr_series |> filter(nRR <= outlier_nRR & x != y)),
       aes(x=nRR,fill=x)) +
  geom_histogram(show.legend = FALSE) +
  labs(x="Normalized RR",y="Count") +
  theme_bw()

df_ordered_nRR <- state_rr_series %>%
  filter(nRR <= outlier_nRR, x != y) %>%
  mutate(x = fct_reorder(x, nRR, .fun = median, .desc = FALSE))  # or FALSE for ascending

ggplot(data=df_ordered_nRR,
       aes(x=nRR,fill=x)) +
  geom_histogram(show.legend = FALSE) +
  labs(x="Normalized RR",y="Count") +
  theme_bw()

ggplot(data = df_ordered_nRR,
       aes(x = nRR, y = x, fill = x)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x = "Normalized RR", y = "State") +
  theme_bw()


ggplot(data=(state_rr_series |> filter(nRR <= outlier_nRR & x != y)),
       aes(x=nRR,y=x,fill=x)) +
  geom_boxplot(show.legend = FALSE) +
  labs(x="Normalized RR",y="State") +
  theme_bw()

aov(nRR ~ x,data = (state_rr_series |> filter(nRR <= outlier_nRR & x != y))) %>% summary()

ggplot(data=(state_rr_series |> filter(x==y)),
       aes(x=date,y=RR,color=x)) +
  geom_line(alpha=0.2) +
  #geom_point(aes(fill=x)) +
  geom_line(y=0,lty=2,size=0.5,color="black") +
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

state_pair_plot("Manitoba", "North Dakota")

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

# Hyper threshold graph code 
# SCALE_FACTOR <- 1
# RR_SIZE <- 1.5
# AXIS_SIZE <- 12

# hmap_state_rr_all <- state_rr_all %>%
#   rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
#   ggplot(aes(x=x,y=y,fill=fill_RR)) +
#   geom_tile() +
#   scale_fill_gradient2(name="RR",
#                        high = "#D67C34",
#                        low = "#4C90C0",
#                        na.value = "black",
#                        limits=c(log10(LB),log10(UB)),
#                        breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
#                        labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB)) +
#   theme_minimal() + 
#   theme(plot.title=element_text(hjust=0.5,size = 16),
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + 
#   theme(axis.text.x = element_text(angle = 45, hjust=1, size = AXIS_SIZE),
#         axis.text.y = element_text(size = AXIS_SIZE)) +
#   labs(title=paste0("All years - Hyper Threshold: ",hCut))


# hmap_state_rr_time <- state_rr_year %>%
#   mutate(x = factor(x,levels=STATE_ORDER)) %>%
#   mutate(y = factor(y,levels=STATE_ORDER)) %>%
#   rowwise %>% mutate(fill_RR = fill_bound(RR)) %>% 
#   ggplot(aes(x=x,y=y,fill=fill_RR)) +
#   facet_wrap(facets=vars(year_date),ncol=5,dir='v') +
#   geom_tile() +
#   scale_fill_gradient2(name="RR",
#                        high = "#D67C34",
#                        low = "#4C90C0",
#                        na.value = "black",
#                        limits=c(log10(LB),log10(UB)),
#                        breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
#                        labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB)) +
#   theme_minimal() + 
#   theme(plot.title=element_text(hjust=0.5,size = 16),
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16)) + 
#   theme(axis.text.x = element_text(angle = 45, hjust=1, size = AXIS_SIZE),
#         axis.text.y = element_text(size = AXIS_SIZE)) +
#   labs(title=paste0("Hyper threshold: ",hCut))

# hmap_state_rr_time