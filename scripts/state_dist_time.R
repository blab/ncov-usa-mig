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
fig_path <- paste0("figs/",SCENARIO,"/time/")

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
                      limits=c(0.0001,0.07),
                      labels=label_number(accuracy = 0.01)) +
  ylim(0,0.5) + 
  labs(title="Normalized Interstate RR over Time",y="nRR",x="Date of RR Snapshot",fill="Median nRR") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(angle=45,hjust=1)
  )
ggsave(filename = paste0(fig_path,"rr_boxplot_time.png"),
       width = 7,
       height = 5,
       units = "in",
       dpi=192)

TEST_DATE <- as.Date("2021-01-01")
test_snap <- state_rr_snap %>% filter(date==TEST_DATE)
ggplot(test_snap %>% filter(x!=y),
       aes(x=euclid_dist,y=nRR)) +
  geom_point() +
  theme_bw()
  
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
       x = "Neighbor Rank (Queen Adjacency)")
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



