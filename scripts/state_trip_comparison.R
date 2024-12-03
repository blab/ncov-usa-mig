#File: state_trip_comparison.R
#Author(s): Amin Bemanian
#Date: 10/2/24
#Description: Comparison of state RRs with interstate travel RRs
#Arguments: 
#--scenario: Scenario corresponding to data files, for this file only use 50 states + DC geography (may use different scenarios for time/strains)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(usmap)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = 'USA', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
state_rr <- fread(fn_rr) %>%
  select(x,y,RR) %>%
  rename(RR_seq = RR)

fn_clust <- paste0("results/",scenario,"/df_state_clusters.tsv")
df_state_cluster <- fread(fn_clust)

df_travel <- data.table::fread("data/interstate_travel_dot.csv") %>%
  rename(x = Origin) %>%
  rename(y = Destination) %>% 
  rename(mode = "Pivot Field Names") %>%
  rename(trips_dir = "Pivot Field Values") %>% #Directional trips
  mutate(trips_dir = as.numeric(trips_dir)) #Switch from into to num since it will confuse the math later

df_travel_air <- data.table::fread("data/interstate_air_travel_dot.csv") %>%
  rename(x = Origin) %>%
  rename(y = Destination) %>% 
  rename(trips_dir = "Pivot Field Values") %>% #Directional trips
  mutate(trips_dir = as.numeric(trips_dir)) #Switch from into to num since it will confuse the math later

STATE_LIST <- unique(df_travel$x)

#Convert from unidirectional to bidirectional trips
df_travel$trips_xy <- 0 
for(i in STATE_LIST){
  for(j in STATE_LIST){
    if(i==j){
      df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir
    }else{
      df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir + df_travel[x == j & y == i]$trips_dir
    }
  }
}

df_travel_air$trips_xy <- 0
for(i in STATE_LIST){
  for(j in STATE_LIST){
    if(i==j){
      df_travel_air[x == i & y == j]$trips_xy <- df_travel_air[x == i & y == j]$trips_dir
    }else{
      df_travel_air[x == i & y == j]$trips_xy <- df_travel_air[x == i & y == j]$trips_dir + df_travel_air[x == j & y == i]$trips_dir
    }
  }
}


df_travel <- df_travel %>%
  group_by(x) %>%
  mutate(trips_x = sum(trips_xy)) %>%
  group_by(y) %>%
  mutate(trips_y = sum(trips_xy)) %>%
  ungroup() %>%
  mutate(trips_total = sum(trips_xy),
         same_state = (x == y), #Use to stratify for out of state travel only
         RR_trips = (trips_xy/trips_x)/(trips_y/trips_total)) %>%
  inner_join(state_rr,by=join_by(x,y))

#Out of state travel only
df_travel_out <- df_travel %>%
  filter(!same_state) %>%
  group_by(x) %>%
  mutate(trips_x = sum(trips_xy)) %>%
  group_by(y) %>%
  mutate(trips_y = sum(trips_xy)) %>%
  ungroup() %>%
  mutate(trips_total = sum(trips_xy),
         RR_trips = (trips_xy/trips_x)/(trips_y/trips_total)) %>%
  rowwise() %>%
  mutate(same_hclust = (df_state_cluster[state==x]$hclust == df_state_cluster[state==y]$hclust))  %>%
  mutate(id_hclust = ifelse(same_hclust, as.character(df_state_cluster[state==x]$hclust),"Different Cluster"))


#Air travel also will be out of state only, also exclude DC
df_travel_air <- df_travel_air %>%
  mutate(same_state = (x == y)) %>%
  filter(!same_state) %>%
  filter(x != "District of Columbia" & y != "District of Columbia") %>%
  inner_join(state_rr,by=join_by(x,y)) %>%
  group_by(x) %>%
  mutate(trips_x = sum(trips_xy)) %>%
  group_by(y) %>%
  mutate(trips_y = sum(trips_xy)) %>%
  ungroup() %>%
  mutate(trips_total = sum(trips_xy),
         RR_trips = (trips_xy/trips_x)/(trips_y/trips_total)) %>%
  rowwise() %>%
  mutate(same_hclust = (df_state_cluster[state==x]$hclust == df_state_cluster[state==y]$hclust)) %>%
  mutate(id_hclust = ifelse(same_hclust, as.character(df_state_cluster[state==x]$hclust),"Different Cluster"))



#Remove double counting for scatter plots
for(i in 1:length(STATE_LIST)){
  df_travel <- filter(df_travel,
                      y != STATE_LIST[i] | x %in% STATE_LIST[1:i])
  df_travel_out <- filter(df_travel_out,
                          y != STATE_LIST[i] | x %in% STATE_LIST[1:i])
  df_travel_air <- filter(df_travel_air,
                          y != STATE_LIST[i] | x %in% STATE_LIST[1:i])
}

rho_travel <- cor(df_travel$RR_seq,df_travel$RR_trips,method = "spearman")
rho_travel.p <- cor.test(df_travel$RR_seq,df_travel$RR_trips,method = "spearman")$p.value

rho_travel_out <- cor(df_travel_out$RR_seq,df_travel_out$RR_trips,method = "spearman")
rho_travel_out.p <- cor.test(df_travel_out$RR_seq,df_travel_out$RR_trips,method = "spearman")$p.value

rho_travel_air <- cor(df_travel_air$RR_seq,df_travel_air$RR_trips,method = "spearman")
rho_travel_air.p <- cor.test(df_travel_air$RR_seq,df_travel_air$RR_trips,method = "spearman")$p.value

plot_travel_seq_all <- ggplot(data = df_travel) +
  geom_point(aes(x=RR_trips,y=RR_seq,color=same_state),alpha=0.4) +
  geom_smooth(aes(x=RR_trips,y=RR_seq),color='firebrick',fill='firebrick',method='gam') +
  scale_x_continuous(transform ='log',
                     name=expression(RR["travel"]),
                     breaks = c(1E-5,1E-3, 1E-1, 1E0, 1E1, 1E3),
                     labels = c(expression(10^{-5}),expression(10^{-3}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{3})),
                     expand = expansion(mult = c(0.18, 0.13))) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.75,1,5,10,20),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.15),10^(1.5))) + 
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom") +
  ggtitle("State Travel versus State Sequences") +
  labs(title="State Travel versus Sequences - All Pairs",
       subtitle = paste0("Spearman's rho: ",round(rho_travel,2)),
       color = "Within Same State")
fn_travel_seq_all <- paste0("figs/",scenario,"/travel_seq_all.jpg")
ggsave(fn_travel_seq_all,plot_travel_seq_all,width = 6,height = 6,dpi = 600)

plot_travel_seq_out <- ggplot(data = df_travel_out) +
  geom_point(aes(x=RR_trips,y=RR_seq,color=same_hclust),alpha=0.4) +
  geom_smooth(aes(x=RR_trips,y=RR_seq),color='firebrick',fill='firebrick',method='gam') +
  scale_x_continuous(transform ='log',
                     name=expression(RR["travel"]),
                     breaks = c(1E-5,1E-3, 1E-1, 1E0, 1E1, 1E3),
                     labels = c(expression(10^{-5}),expression(10^{-3}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{3})),
                     expand = expansion(mult = c(0.18, 0.13))) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.75,1,1.5,2),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.15),10^(0.35))) + 
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom")+
  labs(title="State Travel versus Sequences - Only Different Pairs",
       subtitle = paste0("Spearman's rho: ",round(rho_travel_out,2)),
       color="Within Same Cluster")
fn_travel_seq_out <- paste0("figs/",scenario,"/travel_seq_out.jpg")
ggsave(fn_travel_seq_out,plot_travel_seq_out,width = 6,height=6,dpi = 600)

plot_travel_seq_air <- ggplot(data = df_travel_air) +
  geom_point(aes(x=RR_trips,y=RR_seq,color=same_hclust),alpha=0.4) +
  geom_smooth(aes(x=RR_trips,y=RR_seq),color='firebrick',fill='firebrick',method='gam') +
  scale_x_continuous(transform ='log',
                     name=expression(RR["travel"]),
                     breaks = c(1E-2, 1E-1, 1E0, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13))) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.75,1,1.5,2),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.15),10^(0.35))) + 
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom")+
  labs(title="State Travel versus Sequences - Only Air Travel & Different Pairs",
       subtitle = paste0("Spearman's rho: ",round(rho_travel_air,2)),
       color="Within Same Cluster")
fn_travel_seq_air <- paste0("figs/",scenario,"/travel_seq_air.jpg")
ggsave(fn_travel_seq_air,plot_travel_seq_air,width = 6,height=6,dpi = 600)

#Diagnostic plots, nothing interesting to make into figures
plot_travel_cluster_air <- ggplot(data = df_travel_air) +
  geom_point(aes(x=RR_trips,y=same_hclust),alpha=0.1) +
  scale_x_continuous(transform ='log',
                     name=expression(RR["travel"]),
                     breaks = c(1E-2, 1E-1, 1E0, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)))
  
plot_travel_cluster_out <- ggplot(data = df_travel_out) +
  geom_point(aes(x=RR_trips,y=same_hclust),alpha=0.1) +
  scale_x_continuous(transform ='log',
                     name=expression(RR["travel"]),
                     breaks = c(1E-2, 1E-1, 1E0, 1E1, 1E2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)))

print("Wilcoxon Test for Interstate Data, All Modes - Travel RR vs Hierarchical Cluster ")
print(wilcox.test(RR_trips ~ same_hclust,data=df_travel_out))

print("Wilcoxon Test for Interstate Data, Air Travel Only - Travel RR vs Hierarchical Cluster ")
print(wilcox.test(RR_trips ~ same_hclust,data=df_travel_air))

logistic_trips_hclust_out <- glm(same_hclust ~ log10(RR_trips),family = "binomial",data = filter(df_travel_out, RR_trips > 0))
print(summary(logistic_trips_hclust_out))
logistic_trips_hclust_air <- glm(same_hclust ~ log10(RR_trips),family = "binomial",data = filter(df_travel_air, RR_trips > 0))
print(summary(logistic_trips_hclust_air))

plot_travel_by_cluster_out <- ggplot(data = df_travel_out, aes(x=id_hclust,y=RR_trips)) +
  geom_boxplot() +
  scale_y_continuous(transform='log',
                    name = expression(RR['travel']),
                    breaks = c(1E-3, 1E-2, 1E-1, 1E0, 1E1)) +
  labs(title = "Interstate Travel RR by Cluster ID - All Modes", x = "Cluster ID") +
  theme_classic()
fn_travel_by_cluster_out <- paste0("figs/",scenario,"/travel_cluster_out.jpg")
ggsave(fn_travel_by_cluster_out,plot_travel_by_cluster_out,width = 7,height=5,dpi = 600)

plot_travel_by_cluster_air <- ggplot(data = df_travel_air, aes(x=id_hclust,y=RR_trips)) +
  geom_boxplot() +
  scale_y_continuous(transform='log',
                    name = expression(RR['travel']),
                    breaks = c(1E-3, 1E-2, 1E-1, 1E0, 1E1)) +
  labs(title = "Interstate Travel RR by Cluster ID - Air Travel Only", x = "Cluster ID") +
  theme_classic()
fn_travel_by_cluster_air <- paste0("figs/",scenario,"/travel_cluster_air.jpg")
ggsave(fn_travel_by_cluster_air,plot_travel_by_cluster_air,width = 7,height=5,dpi = 600)

