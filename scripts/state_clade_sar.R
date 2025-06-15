#File: state_clade_sar.R
#Author(s): Amin Bemanian
#Date: 5/5/25
#Description: 
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(tidyverse)
library(data.table)
library(duckdb)
library(dbplyr)
library(sf)
library(spdep)

collect_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
scenario <- "USA"

FREQ_THRESHOLD <- 0.1
DIV_THRESHOLD <- 2
WEEK_WINDOW <- 5 #For SMA smoothing

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)

df_meta <- con %>%
  tbl("metadata")

df_clade <- df_meta %>%
  rename(date_string = date) %>%  # avoid implicit casting
  mutate(
    clean_date_string = sql(
      "
      CASE
        WHEN length(date_string) = 4 AND TRY_CAST(date_string AS INTEGER) IS NOT NULL
        THEN '1900-01-01'
        ELSE date_string
      END
    "
    ),
    safe_date = sql("TRY_CAST(clean_date_string AS DATE)"),
    week = sql("DATE_TRUNC('week', TRY_CAST(clean_date_string AS DATE))")
  ) %>%
  group_by(week, clade_nextstrain, division) %>%
  summarise(n_seq = n(), .groups = "drop") %>%
  arrange(week, division, clade_nextstrain) %>%
  collect() %>%
  filter(week != as.Date("1900-01-01")) #Remove the dummy variables...

zero_grid <- expand_grid(
  week = unique(df_clade$week),
  division = unique(df_clade$division),
  clade_nextstrain = unique(df_clade$clade_nextstrain)
)

df_clade <- left_join(zero_grid,
                      df_clade,
                      by = c("week", "division", "clade_nextstrain")) %>%
  mutate(n_seq = replace_na(n_seq, 0)) %>%
  arrange(week, division, clade_nextstrain) %>%
  group_by(week, division) %>%
  mutate(
    total_seq = sum(n_seq),
    rel_freq = if_else(total_seq > 0, n_seq / total_seq, NA)
  ) %>%
  ungroup %>%
  arrange(week, division, clade_nextstrain)
DBI::dbDisconnect(con) #Everything from here should be local

#Focus time windows during transmission periods 
CLADE_LIST <- unique(df_clade$clade_nextstrain)
DIV_LIST <- unique(df_clade$division)
clades_skipped <- vector()
df_clade_filter <- map_dfr(CLADE_LIST, function(c){
  df_c <- filter(df_clade,clade_nextstrain == c)

  date_range <- df_c %>% 
    filter(rel_freq > FREQ_THRESHOLD) %>% 
    group_by(week) %>%
    summarise(n_div = n_distinct(division), .groups = "drop") %>%
    filter(n_div >= DIV_THRESHOLD)
  
  if (nrow(date_range) == 0) {
    clades_skipped <<- c(clades_skipped, c)
    return(NULL)  # skip this clade
  }
  
  
  date_range <- date_range %>%
    summarise(
      start_date = min(week) - 28,
      end_date = max(week) + 56 #Give more time to watch the decay
    )
  
  df_c %>%
    filter(week >= date_range$start_date & week <= date_range$end_date)
})
message("Skipped clades: ", paste(clades_skipped, collapse = ", "))
sma_inputs <- expand.grid(c = CLADE_LIST, d = DIV_LIST)

#Middle smoothing
#Interval defined as bandwidth/2
smooth_freq<-function(n,t,interval){
  l <- length(n)
  if(l!=length(t)){
    quit("Invalid inputs! n and t must be same length")
  }
  freq <- rep(NA,l)
  for(i in 1:l){
      if(i > interval & i <= (l - interval)){
        n_window <- n[seq(i-interval,i+interval)] %>% replace_na(0) #Replace any NA with 0
        t_window <- t[seq(i-interval,i+interval)] %>% replace_na(0) #Replace any NA with 0
        freq[i] <- sum(n_window)/sum(t_window)
      }
  }
  return(freq)
}

tictoc::tic("Smoothing frequencies")
df_clade_filter <- pmap_dfr(sma_inputs,function(c,d){
  df_d_c <- df_clade_filter %>% 
    filter(division == d) %>%
    filter(clade_nextstrain == c)
#  df_d_c$smooth_freq <- df_d_c$rel_freq %>% 
#    zoo::na.approx(na.rm=FALSE) %>% 
#    zoo::rollmean(k = WEEK_WINDOW,align="center",fill=NA)
  df_d_c$smooth_freq <- smooth_freq(df_d_c$n_seq,
                                    df_d_c$total_seq,
                                    floor(WEEK_WINDOW/2))
  df_d_c$smooth_freq <- replace_na(df_d_c$smooth_freq,0) #Change the leading and tail NAs to 0
  return(df_d_c)
})
tictoc::toc()

write_tsv(df_clade_filter,
          paste0("results/",scenario,"/clade_freq.tsv")
)

#Moran's I calculation
#Load in our shapefile
us_shp <- st_read("data/shp-files/tl_2023_us_state.shp") %>%
  filter(NAME %in% DIV_LIST) %>% #Filters out territories, NAME fortunately is exactly the same as division including for DC
  arrange(NAME)
nb_rook <- poly2nb(us_shp) #Create the neighbor list and weights
lw_rook <- nb2listw(nb_rook, style="W", zero.policy=TRUE)

#Load in a theoretical gravity model attraction matrix, pre-made by calc_grav_matrix.r
df_grav <- fread("data/state_grav.tsv") %>%
  filter(state_x %in% DIV_LIST & state_y %in% DIV_LIST)
grav_wt <- df_grav %>% select(state_x,state_y,grav_wt) %>%
  rename(name = state_x) %>%
  pivot_wider(names_from = state_y, values_from = grav_wt) %>%
  arrange(name) %>%
  select(name, sort(setdiff(names(.), "name"))) %>%
  column_to_rownames("name") %>%
  as.matrix()

inv_dist_wt <- df_grav %>% 
  mutate(inv_dist = 1 / euclid_dist) %>%
  select(state_x,state_y,inv_dist)  %>%
  rename(name = state_x) %>%
  pivot_wider(names_from = state_y, values_from = inv_dist) %>%
  arrange(name) %>%
  select(name, sort(setdiff(names(.), "name"))) %>%
  column_to_rownames("name") %>%
  as.matrix()
  
for(i in 1:length(DIV_LIST)){
  grav_wt[i,i]<-0
  grav_wt[i,] <- grav_wt[i,]/sum(grav_wt[i,])
  inv_dist_wt[i,i]<-0
  inv_dist_wt[i,] <- inv_dist_wt[i,]/sum(inv_dist_wt[i,])
}


#Also load in NHTS data
df_travel <- fread("data/interstate_travel_dot.csv") %>%
  rename(x = Origin) %>%
  rename(y = Destination) %>% 
  rename(mode = "Pivot Field Names") %>%
  rename(trips_dir = "Pivot Field Values") %>% #Directional trips
  mutate(trips_dir = as.numeric(trips_dir))
#Convert from unidirectional to bidirectional trips
df_travel$trips_xy <- 0 
for(i in DIV_LIST){
  for(j in DIV_LIST){
    if(i==j){
      df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir
    }else{
      df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir + df_travel[x == j & y == i]$trips_dir
    }
  }
}

#Calculate a travel RR matrix and convert it to a distance matrix
df_travel <- df_travel %>%
  group_by(x) %>%
  mutate(trips_x = sum(trips_xy)) %>%
  group_by(y) %>%
  mutate(trips_y = sum(trips_xy)) %>%
  ungroup() %>%
  mutate(trips_total = sum(trips_xy),
         same_state = (x == y), #Use to stratify for out of state travel only
         RR = (trips_xy/trips_x)/(trips_y/trips_total),
         trip_wt =(trips_xy)/trips_x)

travel_wt <- df_travel %>% select(c("x","y","trip_wt")) %>%
  rename(name = x) %>%
  pivot_wider(names_from = y, values_from = trip_wt) %>%
  arrange(name) %>%
  select(name, sort(setdiff(names(.), "name"))) %>%
  column_to_rownames("name") %>%
  as.matrix()

for(i in 1:length(DIV_LIST)){
  travel_wt[i,i]<-0
  travel_wt[i,] <- travel_wt[i,]/sum(travel_wt[i,])
}


wt_type <- tibble(t=c("rook","grav","inv_dist","travel"))
moran_inputs <- df_clade_filter %>%
  filter(division=="Washington") %>% #Limit to a single state so you don't have 51 copies
  select(clade_nextstrain,week) %>%
  rename(c = clade_nextstrain) %>%
  rename(w = week) %>% 
  cross_join(wt_type,copy=TRUE) #Add thresholds for using nhts travel weights



#Run a single Moran I test with lots of simulations to get a sense of the residual distribution for the geography
#This will be used to make the non-significant range for the plots
#Selection of clade and week are totally arbitrary except needs to be a week with some non-zero values
#map_bounds <-  df_clade_filter %>% 
#  filter(clade_nextstrain == "21K") %>%
#  filter(week == as.Date("2022-01-03")) %>%
#  left_join(us_shp, ., by=join_by(NAME==division))
#moran_bounds <- moran.mc(x = map_bounds$smooth_freq,
#                         listw = lw_travel,
#                         nsim=1E5,
#                         zero.policy = TRUE,
#                         alternative="two.sided")
#moran_lb <- quantile(moran_bounds$res,0.025)
#moran_ub <- quantile(moran_bounds$res,0.975)
break

df_moran <- pmap_dfr(moran_inputs,function(c,w,t){
  if(t == "rook"){
    lw <- lw_rook
  }else if(t=="grav"){
    lw <- mat2listw(grav_wt,style="W",zero.policy = TRUE)
  }else if(t == "travel"){
    lw <- mat2listw(travel_wt,style="W",zero.policy = TRUE)
  }else if(t == "inv_dist"){
    lw <- mat2listw(inv_dist_wt,style="W",zero.policy = TRUE)
  }else{
    stop("Invalid weight type")
  }
  df_slice <- df_clade_filter %>% 
    filter(clade_nextstrain == c) %>%
    filter(week == w) %>%
    select(division, smooth_freq)
#  map_slice <- left_join(us_shp,
#                         df_slice,
#                        by=join_by(NAME == division))
  mi <- NA
  tryCatch({
    mc_out <- moran.mc(x=df_slice$smooth_freq,
                       listw  = lw,
                       nsim=1, #Skips the simulations
                       zero.policy = TRUE,
                       alternative = "two.sided")
    mi <- mc_out$statistic
  }, error = function(arg){}) #For NA values
  return(tibble(clade_nextstrain = c, 
                week = w,
                wt_type = t,
                moran_i = mi,
                ))},
  .progress = list( #Progress bar for interactive mode
    type = "iterator", 
    format = "Calculating {cli::pb_bar} {cli::pb_percent}",
    clear = TRUE)
)

clade_check <- df_moran %>% filter(clade_nextstrain=="23A")
ggplot(df_moran,
       aes(x = week, y = moran_i, color = wt_type)) +
  geom_line() +
  geom_point() +
  #  geom_ribbon(aes(x=week,ymin=moran_lb,ymax=moran_ub),alpha=0.1) +
  theme_bw() +
  ylim(-1, 1) + 
  facet_wrap(vars(clade_nextstrain), ncol = 6, scales = "free_x") +
  scale_x_date(name = "Calendar Week") +
  scale_y_continuous(name = "Moran's I")

ggsave(filename = paste0("figs/",scenario,"/time/moran_clade.jpg"),
      height = 30,
      width = 30,
      unit = "in",
      dpi = 150)

#Arrival analysis
#Limit to clades that break 50%
target_clades <- df_clade_filter %>% 
  group_by(clade_nextstrain) %>%
  summarize(max_freq = max(smooth_freq),.groups = "drop") %>%
  filter(max_freq > 0.5) %>%
  pull(clade_nextstrain)

arrival_inputs <- expand.grid(c = target_clades, t = wt_type$t)
ARRIVAL_THRESHOLD <- 0.1
df_arrival <- pmap_dfr(arrival_inputs, function(c, t) {
  df_slice <- df_clade_filter %>% filter(clade_nextstrain == c)
  arrival_list <- df_slice %>%
    filter(smooth_freq >= ARRIVAL_THRESHOLD) %>%
    group_by(division) %>%
    summarise(arrival_week = min(week), .groups = "drop") %>%
    right_join(df_slice %>% distinct(division), by = "division") %>%
    arrange(division)
  arrival_list$delta_t <- as.numeric(arrival_list$arrival_week - min(arrival_list$arrival_week, na.rm = TRUE)) / 7
  arrival_list$delta_t[is.na(arrival_list$delta_t)] <- max(arrival_list$delta_t,na.rm=TRUE) + 2
  
  mi <- NA
  p <- NA
  tryCatch({
    valid <- !is.na(arrival_list$delta_t)
    if(t == "rook"){
      x <- arrival_list$delta_t
      lw <- lw_rook
    }else if(t=="grav"){
      x <- arrival_list$delta_t[valid]
      grav_wt_valid <- grav_wt[valid,valid]
      lw <- mat2listw(grav_wt_valid,style="W",zero.policy = TRUE)
    }else if(t == "travel"){
      x <- arrival_list$delta_t[valid]
      travel_wt_valid <- travel_wt[valid,valid]
      lw <- mat2listw(travel_wt_valid,style="W",zero.policy = TRUE)
    }else if(t == "inv_dist"){
      x <- arrival_list$delta_t[valid]
      inv_dist_wt_valid <- inv_dist_wt[valid,valid]
      lw <- mat2listw(inv_dist_wt_valid,style="W",zero.policy = TRUE)
    }else{
      stop("Invalid weight type")
    }
    # Then run moran.mc with na.action = na.exclude
    mc_out <- moran.mc(
      x = x,
      listw = lw,
      nsim = 2e3,
      zero.policy = TRUE,
      na.action = na.exclude,
      alternative = "two.sided"
    )
    mi <- mc_out$statistic
    p <- mc_out$p.value
    lm_out <- localmoran(
      x = x,
      listw = lw,
      zero.policy = TRUE,
      na.action = na.exclude,
      alternative = "two.sided"
    )
  }, error = function(arg) {})
  
  p_signif = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE  ~ " "
  )
  
  return(
    tibble(
      arrival_dates = list(arrival_list),
      local_moran = list(lm_out),
      clade = c,
      wt_type = t,
      moran =  mi,
      moran.p = p,
      moran.signif = p_signif
    )
  )},
  
  .progress = list(#Progress bar for interactive mode
    type = "iterator",
    format = "Calculating {cli::pb_bar} {cli::pb_percent}",
    clear = TRUE)
)

ggplot(df_arrival,aes(x=wt_type,y=clade,fill=moran)) +
  geom_tile(color="white") +
  geom_text(aes(label=moran.signif)) +
  scale_fill_continuous(low="yellow2",high="red3") + theme_bw()

arrival_21C <-left_join(us_shp,df_arrival$arrival_dates[[8]],by=join_by(NAME==division)) 
ggplot(data = arrival_21C %>% filter(!(NAME %in% c("Alaska","Hawaii")))) + 
  ggspatial::geom_sf(aes(fill=delta_t))
