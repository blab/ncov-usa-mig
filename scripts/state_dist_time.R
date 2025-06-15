library(tidyverse)
library(data.table)
library(ggplot2)
library(tictoc)
library(scales)
library(maps)
library(sf)

select <- dplyr::select
df_state_distances <- fread("data/nb_dist_states.tsv")
df_nhts_rr <- fread("data/nhts_rr.csv")

SCENARIO <- "CAM_1000"

attach_distances <- function(df){
  left_join(df,df_state_distances,by=join_by(x==state_x,y==state_y)) %>%
    left_join(df_nhts_rr %>% select(x,y,RR_trips),by=join_by(x==x,y==y)) %>%
    rename(RR_nhts = RR_trips)
}

state_rr_snap <- read_tsv(paste0("results/",SCENARIO,"/time_state/df_state_rr_snap.tsv")) %>%
  attach_distances()
state_rr_series <- read_tsv(paste0("results/",SCENARIO,"/time_state/df_state_rr_series.tsv")) %>%
  attach_distances()

date_medians <- state_rr_snap%>%
  filter(x != y) %>%
  group_by(date) %>%
  summarise(date_median_nRR = median(nRR, na.rm = TRUE)) %>%
  ungroup()

ggplot(state_rr_snap %>%
         filter(x != y) %>%
         left_join(date_medians,by="date"),
       aes(x=date,y=nRR,group=date)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(aes(fill=date_median_nRR),outliers = FALSE,alpha=0.9) +
  scale_fill_gradient(low="gold",high="red2",
                      transform="log10",
                      limits=c(0.001,0.07),
                      labels=label_number(accuracy = 0.01)) +
  ylim(0,0.5) + 
  labs(title="Normalized RR over Time",y="nRR",x="Date of RR Snapshot",fill="Median nRR") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(angle=45,hjust=1)
  )
  
TEST_DATE <- as.Date("2021-01-01")
test_snap <- state_rr_snap %>% filter(date==TEST_DATE)
ggplot(test_snap %>% filter(x!=y),
       aes(x=euclid_dist,y=nRR)) +
  geom_point() +
  theme_bw()
  
fig_path <- paste0("figs/",SCENARIO,"/dist/")
ggplot(state_rr_snap %>% filter(x!=y),
       aes(x=nb_dist,y=nRR,group=interaction(nb_dist,date))) +
  geom_jitter(alpha=0.1) +
  geom_boxplot(aes(fill=as.factor(date)),
               alpha = 0.9,
               outliers = FALSE,
               show.legend = FALSE) +
  theme_bw() +
  ylim(0,0.5) +  
  xlim(0.5,7.5) +
  facet_wrap(vars(date),nrow=2,dir="v") +
  labs(title = "Neighbor Rank and Normalized RR",
       x = "Neighbor Rank (Rook Adjacency)")
ggsave(filename = paste0(fig_path,"nb_boxplot_time.png"),
       width = 2400,
       height = 1200,
       units = "px",
       dpi=192)

ggplot(state_rr_snap %>% filter(x!=y),
       aes(x=euclid_dist,y=nRR)) +
  geom_point(alpha=0.2,
             aes(color=as.factor(date)),
             show.legend = FALSE) +
  geom_smooth(formula = y ~ s(x,k=10,bs="cs")) +
  theme_bw() +
  facet_wrap(vars(date),nrow = 2,dir = "v") +
  xlim(0,5000)  +
  ylim(0,0.5) +
  labs(title = "State Distance and Normalized RR",
       x = "Distance (km)")
ggsave(filename = paste0(fig_path,"euclid_dist_time.png"),
       width = 2400,
       height = 1200,
       units = "px",
       dpi=192)

ggplot(state_rr_snap %>% filter(x!=y),
       aes(x=RR_nhts,y=nRR)) +
  geom_point(alpha=0.2,
             aes(color=as.factor(date)),
             show.legend = FALSE) +
  geom_smooth(formula = y ~ s(x,k=10,bs="cs")) +
  theme_bw() +
  facet_wrap(vars(date),nrow = 2,dir = "v") +
  xlim(0,1)  +
  scale_x_log10() +
  ylim(0,0.5) +
  labs(title = "Travel RR and Normalized RR",
       x = "Travel RR")
ggsave(filename = paste0(fig_path,"nhts_time.png"),
       width = 2400,
       height = 1200,
       units = "px",
       dpi=192)


ggplot(state_rr_snap %>% filter(x!=y),
       aes(x=RR_nhts,y=nRR)) +
  geom_point(alpha=0.2,
             aes(color=as.factor(date)),
             show.legend = FALSE) +
  geom_smooth(formula = y ~ s(x,k=10,bs="cs")) +
  theme_bw() +
  facet_wrap(vars(date),nrow = 2,dir = "v") +
#  scale_y_log10() +
  ylim(0,0.5) +
  scale_x_log10() +
  labs(title = "NHTS Travel RR and Normalized RR",
       x = "RR of Travel")


bump_pairs <- list()
DATE_SNAPS_BUMP <- c("2021-07-01","2022-01-01","2022-07-01","2023-01-01","2023-07-01")
#QUANT_CUTOFF <- 0.7
Z_CUTOFF <- 1

for (i in 1:length(DATE_SNAPS_BUMP)) {
  local_max_set <- state_rr_snap %>% 
    filter(date==DATE_SNAPS_BUMP[i]) %>%
    filter(euclid_dist > 3000 & euclid_dist < 4000)
  mu <- mean(local_max_set$nRR)
  s <- sd(local_max_set$nRR)
  k <- mu + Z_CUTOFF*s
  local_max_set <- filter(local_max_set,nRR >= k)
  bump_pairs[[i]] <- local_max_set %>% select(x,y)  %>%  
    mutate(x_ordered = pmin(x, y), #Remove mirror pairs
           y_ordered = pmax(x, y)) %>%
    distinct(x_ordered, y_ordered, .keep_all = TRUE) %>%
    select(-x_ordered, -y_ordered) 
}

pair_sources <- map_dfr(seq_along(bump_pairs), function(i) {
  bump_pairs[[i]] %>%
    distinct(x, y) %>%          # remove intra-df duplicates
    mutate(source = i)          # tag with source index
})
bump_pair_counts <- pair_sources %>%
  distinct(x, y, source) %>%    # remove repeated source contributions
  count(x, y, name = "n_sources")  # count how many sources

# Get US state map data (includes DC as "district of columbia")
us_states_map <- maps::map("state", fill = TRUE, plot = FALSE)

# Convert to sf and set CRS to WGS84 (EPSG:4326)
us_states_sf <- st_as_sf(us_states_map) %>%
  mutate(ID = tolower(ID)) %>%
  st_set_crs(4326)

# Transform to NAD83 Albers Equal Area (EPSG:5070)
us_states_nad83 <- st_transform(us_states_sf, crs = 5070)

#Map of these bump pairs
state_centroids <- tibble::tibble(
  region = c(state.name,"District of Columbia"),
  lat = c(state.center$y,38.89511),
  lon = c(state.center$x,-77.03637)
) 

bump_pair_coords <- bump_pair_counts %>%
  left_join(state_centroids, by = c("x" = "region")) %>%
  rename(x_lon = lon, x_lat = lat) %>%
  left_join(state_centroids, by = c("y" = "region")) %>%
  rename(y_lon = lon, y_lat = lat) %>%
  filter(!is.na(x_lat) & !is.na(y_lat))
bump_pair_lines <- bump_pair_coords %>%
  rowwise() %>%
  mutate(geometry = st_sfc(
    st_linestring(matrix(c(x_lon, y_lon, x_lat, y_lat), ncol = 2, byrow = FALSE)),
    crs = 4326
  )) %>%
  ungroup() %>%
  st_as_sf() %>%
  st_transform(crs = 5070)
bump_pair_lines$n_sources <- as.numeric(bump_pair_lines$n_sources)

ggplot() +
  geom_sf(data = us_states_nad83, fill = NA, color = "gray10", size = 0.5) +
  geom_sf(data = bump_pair_lines %>% filter(n_sources >= 1), aes(linewidth = n_sources,alpha= n_sources/10), color = "darkred") +
  scale_linewidth_continuous(range = c(0.5, 2.5)) +
  labs(
    title = "Distant Connections",
    linewidth = "Number of Snapshots"
  ) +
  guides(alpha="none")+ 
  theme_void() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

TEST_STATE <- "California"
test_series <- state_rr_series %>% 
  filter(x==TEST_STATE)
ggplot(test_series %>% filter(x != y),
       aes(x=date,y=nRR,color=euclid_dist,group=y)) +
  geom_line(alpha=0.3) + 
  theme_bw() +
  scale_color_continuous(type="viridis",transform="log10") +
  labs(title = "Normalized RR over Time",
       subtitle = paste0("Origin State: ",TEST_STATE),
       x = "Date",
       y = "nRR",
       color = "Neighbor Distance (km)") +
  ylim(0,5)

ggplot(state_rr_series %>% filter(x != y),
       aes(x=date,y=nRR,color=euclid_dist,
           group=interaction(x,y))) +
  geom_line(alpha=0.1) + 
  theme_bw() +
  scale_color_continuous(type="viridis",transform="log10") +
  labs(title = "Normalized RR over Time",
       subtitle = paste0("All states"),
       x = "Date",
       y = "nRR",
       color = "Neighbor Distance (km)") +
  ylim(0,5)

