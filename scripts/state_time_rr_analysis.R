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
REGION_DATA <- fread("data/regions.csv")

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
    select(date, state = x, RR_diag = RR, N_diag = N_pairs) %>%
    mutate(RR_diag = ifelse(N_diag == 0,NA,RR_diag))
  # Join to get RR(x,x) and RR(y,y) for each row
  df_rr <- df_rr %>%
    left_join(rr_diag, by = c("date", "x" = "state")) %>%
    rename(RR_xx = RR_diag) %>%
    left_join(rr_diag, by = c("date", "y" = "state")) %>%
    rename(RR_yy = RR_diag) %>%
    mutate(nRR = RR / rowMeans(across(c(RR_xx, RR_yy)), na.rm = FALSE))
return(df_rr)
}

fn_db <- paste0("db_files/db_",SCENARIO,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db,read_only = TRUE)
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

#Snapshot analysis will be every quarter (3 months)
#Generalize upper and lowerbound definition
TIME_LB_YEARS <- seq(as.Date("2020-01-01"),as.Date("2025-01-01"),by="3 mo")
TIME_UB_YEARS <- TIME_LB_YEARS[2:length(TIME_LB_YEARS)]-1 #Drop the first bound
TIME_LB_YEARS <- TIME_LB_YEARS[1:length(TIME_LB_YEARS)-1] #Drop the last bound

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

MONTH_BUFFER = 3 
# Make date sequences for the time series analysis
mid_start_date <- as.Date("2020-03-01")
mid_end_date <- as.Date("2024-10-01")
mid_date_vector <- seq(from = mid_start_date, to = mid_end_date,by = "4 weeks")
lb_date_vector <- mid_date_vector - 28 * MONTH_BUFFER
ub_date_vector <- mid_date_vector + 28 * MONTH_BUFFER


tic("Series anaylsis")
state_rr_series <- NULL
for(i in 1:length(lb_date_vector)){
  if(i == 1){
    state_rr_series <- con %>% 
      bind_pairs_exp("division",time_bounds = c(lb_date_vector[i],ub_date_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(lb_date_vector[i],ub_date_vector[i])))
  } else{
    rr_series <- con %>% 
      bind_pairs_exp("division",time_bounds = c(lb_date_vector[i],ub_date_vector[i])) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(date = mean(c(lb_date_vector[i],ub_date_vector[i])))
    state_rr_series <- bind_rows(state_rr_series,rr_series)
  } 
}
toc()

state_region_list <- REGION_DATA %>%
  select(state,bea_reg,country)

state_rr_series <- state_rr_series %>%
  left_join(state_region_list, join_by(x == state)) %>%
  rename(region_x = bea_reg) %>%
  rename(country_x = country) %>%
  left_join(state_region_list, join_by(y == state)) %>%
  rename(region_y = bea_reg) %>%
  rename(country_y = country) %>%
  mutate(
    pair_type = case_when(
      country_x != country_y        ~ "International",
      region_x  != region_y         ~ "Different Regions, Same Country",
      x != y                        ~ "Different States, Same Region",
      TRUE                          ~ "Same State/Province"
  )) %>%
  mutate(pair_type = factor(pair_type,
                            levels = c("Same State/Province",
                                       "Different States, Same Region",
                                       "Different Regions, Same Country",
                                       "International")))

state_rr_series <- normalized_state_rr(state_rr_series)
write_tsv(state_rr_series,file=paste0(fn_out_path,"df_state_rr_series.tsv"))
dbDisconnect(con)


PAIR_LIST <- data.frame(
  x = c("Mexico",
        "Mexico",
        "Mexico",
        "Mexico",
        "Ontario",
        "Ontario",
        "Ontario",
        "Ontario",
        "California",
        "California",
        "California",
        "California"),
  y = c("Texas",
        "California",
        "Illinois",
        "New York",
        "New York",
        "Quebec",
        "British Columbia",
        "Mexico",
        "New York",
        "Illinois",
        "Texas",
        "Washington")
)

sub_RR_snap <- state_rr_snap %>%
  inner_join(PAIR_LIST, by=c("x","y")) %>%
  mutate(pair_name = paste0(x,"-",y)) %>%
  arrange(pair_name)

ggplot(sub_RR_snap,
       aes(x=date,y=pair_name,fill=nRR)) +
  geom_tile(color="#555") +
  scale_x_date(date_breaks = "6 months") +
  labs(x = "Date", y = "State/Province Pair") +
  scale_fill_gradient(low="#EEE",
                      high="#F73") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size= 11)) +
    theme_minimal()


unique_series <- state_rr_series %>% 
  mutate(duplicate = x >= y) %>% 
  filter (duplicate != TRUE) %>%
  mutate(pair_name = paste0(x,"-",y)) %>%
  arrange(pair_name)

ggplot(unique_series,
       aes(x = date, y = nRR, group = date, colour = pair_type)) +
  geom_boxplot(alpha = 0.2) +
  geom_smooth(aes(group = 1)) +
  coord_cartesian(ylim=c(0,0.5)) +
    facet_wrap(~ pair_type, ncol = 1) +
  theme_bw() +
  theme(legend.position = "none")

ts_nRR_df <- unique_series %>%
  mutate(
    date_num  = as.numeric(date),
    pair_type = factor(pair_type),
    pair_name = factor(pair_name)
  ) %>%
  filter(!is.na(nRR), !is.na(date_num), !is.na(pair_type), !is.na(pair_name))


ts_nRR_model <- mgcv::gam(
  nRR ~ s(date_num,bs="cs",k=10) + pair_type + s(pair_name,bs="re"),
  data=ts_nRR_df
)

summary(ts_nRR_model)
mgcv::gam.check(ts_nRR_model)
AIC(ts_nRR_model)
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

state_pair_plot("New York", "California")

