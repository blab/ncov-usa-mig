#File: state_trip_comparison.R
#Author(s): Amin Bemanian
#Date: 10/2/24
#Description: Comparison of state RRs with interstate travel RRs
#Arguments: 
#--scenario: Scenario corresponding to data files, for this file only use 50 states + DC geography (may use different scenarios for time/strains)

library(argparse)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggrepel)
library(usmap)
library(scales)
library(rsample)    # vfold_cv()
library(yardstick)  # rmse()
library(broom)      # tidy / glance

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = 'USA', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
scenario <- "CAM_1000"

fn_rr <- paste("results/",scenario,"/df_RR_by_state.tsv",sep="")
state_rr <- fread(fn_rr) %>%
  select(x,y,euclid_dist,min_cbsa_dist,RR) %>%
  rename(RR_seq = RR)

df_travel <- data.table::fread("data/interstate_travel_dot.csv") %>%
  rename(x = Origin) %>%
  rename(y = Destination) %>% 
  rename(mode = "Pivot Field Names") %>%
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
  rowwise()

df_travel_air <- read_csv("data/rr_air_db1b.csv") %>%
  inner_join(state_rr,by=join_by(x,y))

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

rho_travel_air <- cor(df_travel_air$RR_seq,df_travel_air$RR_air,method = "spearman")
rho_travel_air.p <- cor.test(df_travel_air$RR_seq,df_travel_air$RR_air,method = "spearman")$p.value

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
                     breaks = c(0.01,0.1,0.5,1,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     labels  = label_number(accuracy = 0.01),
                     limits = c(10^(-1.5),10^(2.5))) + 
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom") +
  labs(title="NHTS Travel versus Identical Sequences - All Pairs",
       subtitle = paste0("Spearman's rho: ",round(rho_travel,2)),
       color = "Within Same State")
fn_travel_seq_all <- paste0("figs/",scenario,"/travel_seq_all.jpg")
ggsave(fn_travel_seq_all,plot_travel_seq_all,width = 6,height = 4,dpi = 600)

plot_travel_seq_out <- ggplot(data = df_travel_out) +
  geom_point(aes(x=RR_trips,y=RR_seq),alpha=0.4) +
  geom_smooth(aes(x=RR_trips,y=RR_seq),color='firebrick',fill='firebrick',method='gam') +
  scale_x_continuous(transform ='log',
                     name=expression(RR["travel"]),
                     breaks = c(1E-5,1E-3, 1E-1, 1E0, 1E1, 1E3),
                     labels = c(expression(10^{-5}),expression(10^{-3}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{3})),
                     expand = expansion(mult = c(0.18, 0.13))) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.1,0.5,1,2,10),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-1.1),10^(1.1))) + 
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom")+
  labs(title="NHTS Travel versus Identical Sequences - Only Different States",
       subtitle = paste0("Spearman's rho: ",round(rho_travel_out,2)))
fn_travel_seq_out <- paste0("figs/",scenario,"/travel_seq_out.jpg")
ggsave(fn_travel_seq_out,plot_travel_seq_out,width = 6,height=4,dpi = 600)

plot_travel_seq_air <- ggplot(data = df_travel_air) +
  geom_point(aes(x=RR_air,y=RR_seq),alpha=0.4) +
  geom_smooth(aes(x=RR_air,y=RR_seq),color='firebrick',fill='firebrick',method='gam') +
  scale_x_continuous(transform ='log',
                     name=expression(RR["air travel"]),
                     breaks = c(1E-5,1E-3, 1E-1, 1E0, 1E1, 1E3),
                     labels = c(expression(10^{-5}),expression(10^{-3}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{3})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-5),10^(3))) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.01,0.1,0.5,1,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     labels  = label_number(accuracy = 0.01),
                     limits = c(10^(-1.5),10^(2.5))) + 
  geom_hline(yintercept = 1) + #Reference point
  theme_classic() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom")+
  labs(title="Air Travel versus Sequences",
       subtitle = paste0("Spearman's rho: ",round(rho_travel_air,2)))
fn_travel_seq_air <- paste0("figs/",scenario,"/travel_seq_air.jpg")
ggsave(fn_travel_seq_air,plot_travel_seq_air,width = 6,height=4,dpi = 600)


## Cross validation work by ChatGPT to determine thresholding for air
cv_score <- function(df, short_cut, long_cut,
                     three_level = TRUE,            # FALSE = 2-level
                     K = 10,                        # k-fold
                     metric = rmse) {               # any yardstick fn
  # 1. slice the data -------------------------------------------------
  if(three_level){
    df2 <- df %>%
      mutate(
        flight_length = case_when(
          min_cbsa_dist  <  short_cut ~ "Short",
          min_cbsa_dist  > long_cut   ~ "Long",
          TRUE ~ "Medium"
        )
      ) %>%
      mutate(flight_length = factor(flight_length,levels=c("Short","Medium","Long")))
  }else {
    df2 <- df %>%
      mutate(flight_length = case_when(
        min_cbsa_dist<short_cut ~ "Short",
        TRUE ~ "Long")) %>%
      mutate(flight_length = factor(flight_length,levels=c("Short","Long")))
  }
  
  # 2. build folds ----------------------------------------------------
  folds <- vfold_cv(df2, v = K, strata = flight_length)
  
  # 3. fit & score each fold -----------------------------------------
  fold_stats <- map_dfr(folds$splits, function(split) {
    train <- analysis(split)
    test  <- assessment(split)
    mod  <- lm(log10(RR_seq) ~ log10(RR_air) : flight_length, data = train)
    pred <- tibble(
      truth    = log10(test$RR_seq),          # now a real column
      estimate = predict(mod, test)
    )
    
    tibble(rmse = rmse(pred, truth, estimate)$.estimate)
  })
  
  mean_rmse <- mean(fold_stats$rmse)
  print(paste0("Short cut: ",short_cut))
  if(three_level){print(paste0("Long cut: ",long_cut))}
  print(paste0("Mean RMSE: ",round(mean_rmse,3)))
  return(mean_rmse)
}

# --------------------------------------------------------------------
# Build the distance grid to search ----------------------------------
grid <- expand.grid(
  short_cut = seq(100, 1000, by = 25),
  long_cut  = seq(200, 5000, by = 25)
) %>%
  filter(short_cut < long_cut)

# --------------------------------------------------------------------
# 1) Three-level optimisation ----------------------------------------
three_results <- grid %>%
  mutate(rmse = pmap_dbl(
    list(short_cut, long_cut),
    cv_score,
    df           = df_travel_air %>% filter(x != y),
    three_level  = TRUE,
    K            = 10
  ))

best3 <- slice_min(three_results, rmse, n = 1)
best3

# --------------------------------------------------------------------
# 2) Two-level optimisation ------------------------------------------
short_cut <- seq(100,5000, by = 25)
two_results <- grid %>%
  mutate(rmse = pmap_dbl(
    list(short_cut, long_cut),
    cv_score,
    df  = df_travel_air %>% filter(x != y),
    three_level  = FALSE,
    K            = 10
  ))

best2 <- slice_min(two_results, rmse, n = 1)
best2
#> short_cut long_cut   rmse
#>      150      600   0.135

# --------------------------------------------------------------------
# Optional: final refits with the best thresholds --------------------
make_fitted_df <- function(df, short_cut, long_cut, three_level = TRUE) {
  df2 <- df %>%
    mutate(
      flight_length = case_when(
        min_cbsa_dist <  short_cut ~ "Short",
        euclid_dist    > long_cut  ~ "Long",
        TRUE                        ~ "Medium"
      )
    )
  
  if (!three_level) {
    df2 <- df2 %>% filter(flight_length != "Medium")
  }
  df2$flight_length <- droplevels(df2$flight_length)
  df2
}

df3 <- make_fitted_df(df_travel_air %>% filter(x != y) %>%
                        left_join(df_travel %>% select(x, y, RR_trips),
                                  join_by(x, y)),
                      best3$short_cut, best3$long_cut, TRUE)

df2 <- make_fitted_df(df_travel_air %>% filter(x != y) %>%
                        left_join(df_travel %>% select(x, y, RR_trips),
                                  join_by(x, y)),
                      best2$short_cut, best2$long_cut, FALSE)

model3 <- lm(log10(RR_seq) ~ log10(RR_air) : flight_length, data = df3)
model2 <- lm(log10(RR_seq) ~ log10(RR_air) : flight_length, data = df2)

glance(model3)$adj.r.squared
glance(model2)$adj.r.squared

df_travel_dist <- df_travel_air %>%
  filter(x != y) %>%
  left_join(df_travel %>% select(x,y,RR_trips),join_by(x,y)) %>%
  mutate(flight_length = ifelse(min_cbsa_dist < 200, "Short",
                            ifelse(euclid_dist > 500, "Long","Medium"))) %>%
  mutate(flight_length = factor(flight_length,levels = c("Short","Medium","Long")))


travel_comps <- list(
  c("Short",  "Medium"), 
  c("Medium", "Long")
)


ggplot(df_travel_dist, aes(x = flight_length, y = RR_seq)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(                       # <-- adds the bars
    comparisons = travel_comps,                 # list of 2-element vectors
    label      = "p.signif",                # show “*, **, ns”, etc.
    method     = "wilcox.test"              # or "t.test", "kruskal.test", …
  ) +
  theme_bw() +
  labs(x = "Travel Length") +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.01,0.1,0.5,1,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     labels  = label_number(accuracy = 0.01),
                     limits = c(10^(-1.5),10^(1.5))) + 
  geom_hline(yintercept = 1) 

ggplot(df_travel_dist, aes(x = flight_length, y = RR_trips)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(                       # <-- adds the bars
    comparisons = travel_comps,                 # list of 2-element vectors
    label      = "p.signif",                # show “*, **, ns”, etc.
    method     = "wilcox.test"              # or "t.test", "kruskal.test", …
  ) +
  theme_bw() +
  labs(x = "Travel Length") +
  scale_y_continuous(transform ='log',
                     name=expression(RR["NHTS"]),
                     breaks = c(0.01,0.1,0.5,1,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     labels  = label_number(accuracy = 0.01),
                     limits = c(10^(-1.5),10^(1.5))) + 
  geom_hline(yintercept = 1) 
  
ggplot(df_travel_dist) +
  geom_point(aes(x=RR_trips,y=RR_seq),color="navy",alpha = 0.2)  +
  geom_smooth(aes(x=RR_trips,y=RR_seq),color="navy",alpha = 0.8) +
  scale_x_continuous(transform ='log',
                     name=expression(RR["air travel"]),
                     breaks = c(1E-5,1E-3, 1E-1, 1E0, 1E1, 1E3),
                     labels = c(expression(10^{-5}),expression(10^{-3}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{3})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-5),10^(3))) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.01,0.1,0.5,1,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     labels  = label_number(accuracy = 0.01),
                     limits = c(10^(-1.5),10^(2.5))) + 
  geom_point(aes(x=RR_air,y=RR_seq),color="firebrick",alpha = 0.2)  +
  geom_smooth(aes(x=RR_air,y=RR_seq),color="firebrick",alpha = 0.8) +
  scale_x_continuous(transform ='log',
                     name=expression(RR["travel"]),
                     breaks = c(1E-5,1E-3, 1E-1, 1E0, 1E1, 1E3),
                     labels = c(expression(10^{-5}),expression(10^{-3}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{3})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-5),10^(3))) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(0.01,0.1,0.5,1,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     labels  = label_number(accuracy = 0.01),
                     limits = c(10^(-1.5),10^(2.5))) + 
  geom_hline(yintercept = 1) + #Reference point
  theme_bw() + 
  theme(plot.title=element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),

        axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "bottom")+
  facet_wrap(vars(flight_length),nrow = 3)


zero_correction <- min(df_travel_dist[df_travel_dist$RR_trips > 0,'RR_trips'])

lm(data = df_travel_dist,
   formula = log10(RR_seq) ~ log10(RR_trips+zero_correction)) %>% AIC()

lm(data = df_travel_dist,
   formula = log10(RR_seq) ~ log10(RR_air)) %>% AIC()

lm(data = df_travel_dist,
   formula = log10(RR_seq) ~ log10(RR_air) : flight_length) %>% AIC()

lm(data = df_travel_dist,
   formula = log10(RR_seq) ~ euclid_dist) %>% AIC()

lm(data = df_travel_dist,
   formula = log10(RR_seq) ~ min_cbsa_dist) %>% AIC()

lm(data = df_travel_dist,
   formula = log10(RR_seq) ~ log10(RR_trips+zero_correction) + log10(RR_air) : flight_length) %>% summary()
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

