#File: state_regression.R
#Author(s): Amin Bemanian
#Date: 10/2/24
#Description: Comparison of state RRs with interstate travel RRs and regression formation
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
library(sandwich)
library(lmtest)
library(splines)


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = 'CAM_1000', help = 'Which scenario to perform the analysis on')
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


df_move <- read_tsv("data/safegraph_states_adj_pullano.tsv") %>%
  rename(x = state_origin,
         y = state_destination,
         n_move_xy = n_movement_origin_destination,
         ) %>%
  mutate(n_move_all = sum(n_move_xy,na.rm = TRUE)) %>% 
  group_by(x) %>%
  mutate(n_move_x = sum(n_move_xy,na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(y) %>%
  mutate(n_move_y = sum(n_move_xy,na.rm = TRUE)) %>%
  ungroup()

df_move <- df_move %>% 
  left_join(
    df_move %>%
      select(x_rev = x, y_rev = y, n_move_xy) %>%
      rename(n_move_yx = n_move_xy),
    by = c("x" = "y_rev", "y" = "x_rev")
  ) %>%
  mutate(n_move_avg = 1/2 * (n_move_xy+n_move_yx)) %>%
  mutate(RR_move = n_move_avg * n_move_all/(n_move_x * n_move_y)) %>%
  full_join(df_travel,by=join_by(x,y)) #Join with travel so we have all the values

ggplot(df_move) + 
  aes(x=RR_trips,y=RR_move) + 
  geom_point(alpha = 0.4,aes(color=same_state)) + 
  geom_abline(slope=1,color="red",linetype="dashed") +
  scale_x_continuous(transform = "log10",name = expression(RR["NHTS"]),limits=c(1E-3,1E3)) + 
  scale_y_continuous(transform="log10",name = expression(RR["Safegraph"]),limits=c(1E-3,1E3)) +
  theme_bw()

num_states <- nrow(df_move) %>% sqrt()
k <- 1:(num_states^2)
i <- (k-1) %% (num_states) + 1
j <- floor(k/num_states) + 1
k_sym <- num_states * (i-1) + j

#Remove double counting for scatter plots
for(i in 1:length(STATE_LIST)){
  df_travel <- filter(df_travel,
                      y != STATE_LIST[i] | x %in% STATE_LIST[1:i])
  df_travel_out <- filter(df_travel_out,
                          y != STATE_LIST[i] | x %in% STATE_LIST[1:i])
  df_travel_air <- filter(df_travel_air,
                          y != STATE_LIST[i] | x %in% STATE_LIST[1:i])
  df_move <- filter(df_move,
                          y != STATE_LIST[i] | x %in% STATE_LIST[1:i])
}

rho_travel <- cor(df_travel$RR_seq,df_travel$RR_trips,method = "spearman")
rho_travel.p <- cor.test(df_travel$RR_seq,df_travel$RR_trips,method = "spearman")$p.value

rho_travel_out <- cor(df_travel_out$RR_seq,df_travel_out$RR_trips,method = "spearman")
rho_travel_out.p <- cor.test(df_travel_out$RR_seq,df_travel_out$RR_trips,method = "spearman")$p.value

rho_travel_air <- cor(df_travel_air$RR_seq,df_travel_air$RR_air,method = "spearman")
rho_travel_air.p <- cor.test(df_travel_air$RR_seq,df_travel_air$RR_air,method = "spearman")$p.value

rho_move <- cor(df_move$RR_seq,df_move$RR_move,method = "spearman",use = "complete.obs")
rho_move.p <- cor.test(df_move$RR_seq,df_move$RR_move,method = "spearman",use = "complete.obs")$p.value

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
  labs(title="NHTS Travel versus Identical Sequences",
       subtitle = paste0("Spearman's rho: ",round(rho_travel,2)),
       color = "Within Same State")
fn_travel_seq_all <- paste0("figs/",scenario,"/travel_seq_all.jpg")
ggsave(fn_travel_seq_all,plot_travel_seq_all,width = 6,height = 4,dpi = 192)

plot_move_seq <- ggplot(data = df_move) +
  geom_point(aes(x=RR_move,y=RR_seq,color=same_state),alpha=0.4) +
  geom_smooth(aes(x=RR_move,y=RR_seq),color='firebrick',fill='firebrick',method='gam') +
  scale_x_continuous(transform ='log',
                     name=expression(RR["mobility"]),
                     breaks = c(1E-5,1E-3,1E-1,1E0, 1E1, 1E3),
                     labels = c(expression(10^-5),expression(10^{-3}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{3})),
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
  labs(title="Mobility Score versus Identical Sequences",
       subtitle = paste0("Spearman's rho: ",round(rho_move,2)),
       color = "Within Same State")
fn_move_seq <- paste0("figs/",scenario,"/move_seq.jpg")
ggsave(fn_move_seq,plot_move_seq,width = 6,height = 4,dpi = 192)

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
ggsave(fn_travel_seq_out,plot_travel_seq_out,width = 6,height=4,dpi = 192)

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
ggsave(fn_travel_seq_air,plot_travel_seq_air,width = 6,height=4,dpi = 192)

OPTIMIZE_AIR <- FALSE
if(OPTIMIZE_AIR){
  ## Cross validation work by ChatGPT to determine thresholding for air
  cv_score <- function(df, short_cut, long_cut,
                       three_level = TRUE,            # FALSE = 2-level
                       K = 10,                        # k-fold
                       metric = rmse) {               # any yardstick fn
    # 1. slice the data -------------------------------------------------
    if(three_level){
      df_test <- df %>%
        mutate(
          flight_length = case_when(
            min_cbsa_dist  <  short_cut ~ "Short",
            min_cbsa_dist  > long_cut   ~ "Long",
            TRUE ~ "Medium"
          )
        ) %>%
        mutate(flight_length = factor(flight_length,levels=c("Short","Medium","Long")))
    }else {
      df_test <- df %>%
        mutate(flight_length = case_when(
          min_cbsa_dist<short_cut ~ "Short",
          TRUE ~ "Long")) %>%
        mutate(flight_length = factor(flight_length,levels=c("Short","Long")))
    }
    
    # 2. build folds ----------------------------------------------------
    folds <- vfold_cv(df_test, v = K, strata = flight_length)
    
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
    short_cut = seq(100, 2000, by = 25),
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
  ggplot(three_results, aes(x = short_cut, y = long_cut,fill=rmse)) + 
    geom_tile() + 
    scale_fill_viridis_c() +
    labs(x = "Short Threshold",
           y = "Long Threshold") +
    theme_bw()
  best3 <- slice_min(three_results, rmse, n = 1)
  best3
  
  # --------------------------------------------------------------------
  # 2) Two-level optimisation ------------------------------------------
  grid<- expand.grid(
    short_cut = seq(50,5000, by = 10),
    long_cut = 5000
  )
  two_results <- grid %>%
    mutate(rmse = pmap_dbl(
      list(short_cut, long_cut),
      cv_score,
      df  = df_travel_air %>% filter(x != y),
      three_level  = FALSE,
      K            = 10
    ))
  ggplot(two_results, aes(x = short_cut, y = rmse)) + 
    geom_line(aes(color=rmse)) +
    geom_point(aes(color=rmse)) + 
    
    scale_color_viridis_c() +
    labs(x = "Distance Threshold") +
    theme_bw()
    
  best2 <- slice_min(two_results, rmse, n = 1)
  best2
}

#Set thresholds based on optimization, hardcoding so it does not have to be run each time
df_travel_dist3 <- df_travel_air %>%
  filter(x != y) %>%
  left_join(df_travel %>% select(x,y,RR_trips),join_by(x,y)) %>%
  left_join(df_move %>% select(x,y,RR_move),join_by(x,y)) %>%
  mutate(flight_length = ifelse(min_cbsa_dist < 300, "Short",
                            ifelse(min_cbsa_dist > 1800, "Long","Medium"))) %>%
  mutate(flight_length = factor(flight_length,levels = c("Short","Medium","Long")))

#Single threshold version
df_travel_dist2 <- df_travel_air %>%
  filter(x != y) %>%
  left_join(df_travel %>% select(x,y,RR_trips),join_by(x,y)) %>%
  left_join(df_move %>% select(x,y,RR_move),join_by(x,y)) %>%
  mutate(flight_length = ifelse(min_cbsa_dist < 300, "Short",
                                "Long")) %>%
  mutate(flight_length = factor(flight_length,levels = c("Short","Long")))


df_travel_dist <- df_travel_dist3

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
  
epsilon <- 1E-4

#Make a custom df with rho by flight length
df_air_plot <- df_travel_dist %>%
  group_by(flight_length) %>%
  mutate(
    rho = cor(RR_air + epsilon, RR_seq + epsilon,
              method = "spearman", use = "complete.obs"),
    facet_label = paste0(flight_length, " (rho = ", round(rho, 2), ")")
  ) %>%
  ungroup() %>%
  mutate(facet_label = factor(
    facet_label,
    levels = c(
      paste0("Short (rho = ", round(unique(rho[flight_length == "Short"]), 2), ")"),
      paste0("Medium (rho = ", round(unique(rho[flight_length == "Medium"]), 2), ")"),
      paste0("Long (rho = ", round(unique(rho[flight_length == "Long"]), 2), ")")
    )
  ))

plot_air_length <- ggplot(df_air_plot) +
  geom_point(aes(x = RR_air + epsilon, y = RR_seq + epsilon),
             color = "black", alpha = 0.2) +
  geom_smooth(aes(x = RR_air + epsilon, y = RR_seq + epsilon),
              color = "firebrick", method = "gam", alpha = 0.8) +
  scale_x_continuous(transform = 'log',
                     name = expression(RR["air travel"]),
                     breaks = 10^(-2:2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^-2, 10^1)) +
  scale_y_continuous(transform = 'log',
                     name = expression(RR["identical sequences"]),
                     breaks = 10^(-2:2),
                     labels = c(expression(10^{-2}),expression(10^{-1}),expression(10^{0}),expression(10^{1}),expression(10^{2})),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^-1.2, 10^1.2)) +
  geom_hline(yintercept = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  facet_wrap(vars(facet_label), ncol = 3)
ggsave(plot_air_length,
       filename = paste0("figs/",scenario,"/travel_seq_air_length.jpg"),
       width = 15,
       height = 6,
       dpi = 192
)

zero_correction <- min(df_travel_dist[df_travel_dist$RR_trips > 0,'RR_trips'])

#Modified summary to add specific parameters of interest
model_summary <- function(m,calc_bootstrap=TRUE) {
  m %>% summary() %>% print()
  
  # Count the number of predictors
  n_terms <- length(attr(terms(m), "term.labels"))
  
  if (n_terms >= 2) {
    cat("\nVIF:\n")
    print(car::vif(m))
  } else {
    cat("\nVIF not shown: model has fewer than 2 predictors.\n")
  }
  
  cat("\n Breusch-Pagan Test for Heteroskedasticity: \n")
  bptest(m) %>% print()
  plot(m,which=1)
  
  # Robust coefficient test
  cat("\nRobust Coefficient Test (HC1 robust SEs):\n")
  coeftest(m, vcov = sandwich::vcovHC(m, type = "HC1")) %>% print()
  
  if (calc_bootstrap) {
    cat("\nBootstrapped CIs for Coefficients:\n")
    boot_m <- car::Boot(m, R = 5000)
    print(summary(boot_m))
    print(confint(boot_m))
  }
  cat("AIC:\n")
  m %>% AIC() %>% cat("\n")
}

#Function for observed vs predicted plots for RR
#m - Model
#RR_obs - Vector of the observed RR
#model_name - To be used for title
#log_predict - Whether the output is log-transformed
#epsilon for log-safety
plot_pred_RR <- function(m, RR_obs, model_name, log_predict = TRUE, eps = 1e-10) {
  
  yhat <- predict(m)
  if(log_predict){yhat <- 10^yhat} #Correction for log-transformation
  
  df <- tibble(RR_obs = as.numeric(RR_obs),
               RR_hat = as.numeric(yhat)) %>%
    tidyr::drop_na() %>% 
    mutate(RR_obs = RR_obs + eps,
           RR_hat = RR_hat + eps)
  
  
  rmse_val <- yardstick::rmse_vec(truth = df$RR_obs, estimate = df$RR_hat)
  print(rmse_val)
  rmse_val <- as.numeric(rmse_val)
  range_vals <- range(c(df$RR_obs, df$RR_hat), na.rm = TRUE)
  
  ggplot(df, aes(x = RR_hat, y = RR_obs)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    coord_equal() +
    scale_x_continuous(transform = "log10", name = "Predicted") +
    scale_y_continuous(transform = "log10", name = "Observed") +
    ggtitle(sprintf("%s - RMSE: %.3f", model_name, rmse_val)) +
    coord_equal(xlim = range_vals, ylim = range_vals) +
    theme_bw()
}


df_model_data <- df_travel_dist %>% filter(x != y & !is.na(RR_move))
y_obs <- df_model_data$RR_seq

lm_trips <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ log10(RR_trips+zero_correction)) 
lm_trips %>% model_summary()
lm_trips_plot <- plot_pred_RR(lm_trips,y_obs,"NHTS Trips")

slm_trips <- lm(data = df_model_data,
     formula = log10(RR_seq) ~ bs(log10(RR_trips+zero_correction),df=6)) 
slm_trips %>% model_summary()
slm_trips_plot <- plot_pred_RR(slm_trips,y_obs,"NHTS Trips - Spline")

lm_move <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ log10(RR_move))
lm_move %>% model_summary()
lm_move_plot <- plot_pred_RR(lm_move,y_obs,"Safegraph Mobility")

slm_move <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ bs(log10(RR_move),df=6))
slm_move %>% model_summary(calc_bootstrap = FALSE)
slm_move_plot <- plot_pred_RR(slm_move,y_obs,"Safegraph Mobility - Spline")

lm_air <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ log10(RR_air)) 
lm_air %>% model_summary()
lm_air_plot <- plot_pred_RR(lm_air,y_obs,"Air Travel")

slm_air <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ bs(log10(RR_air),df=6))
slm_air %>% model_summary(calc_bootstrap = FALSE)
slm_air_plot <- plot_pred_RR(slm_air,y_obs,"Air Travel - Spline")

lm_air_int <- lm(data = df_model_data,
                 formula = log10(RR_seq) ~ log10(RR_air) : flight_length)
lm_air_int %>% model_summary()
lm_air_int_plot <- plot_pred_RR(lm_air_int,y_obs,"Air Travel : Flight Length")


slm_air_int <- lm(data = df_model_data,
                  log10(RR_seq) ~ bs(log10(RR_air), df = 6) * bs(min_cbsa_dist,df=6))
slm_air_int %>% model_summary(calc_bootstrap = FALSE)
slm_air_int_plot <- plot_pred_RR(slm_air_int,y_obs,"Air Travel x CBSA Distance - Spline")

lm_euclid <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ euclid_dist)
lm_euclid %>% model_summary()
lm_euclid_plot <- plot_pred_RR(lm_euclid,y_obs,"Euclidean Distance")

slm_euclid <- lm(data = df_model_data,
                formula = log10(RR_seq) ~ bs(euclid_dist,df=6))
slm_euclid %>% model_summary(calc_bootstrap = FALSE)
slm_euclid_plot <- plot_pred_RR(slm_euclid,y_obs,"Euclidean Distance - Spline")

lm_cbsa <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ min_cbsa_dist)
lm_cbsa %>% model_summary()
lm_cbsa_plot <- plot_pred_RR(lm_cbsa,y_obs,"CBSA Distance")

slm_cbsa <- lm(data = df_model_data,
              formula = log10(RR_seq) ~ bs(min_cbsa_dist,df=6))
slm_cbsa %>% model_summary(calc_bootstrap = FALSE)
slm_cbsa_plot <- plot_pred_RR(slm_cbsa,y_obs,"CBSA Distance")

#Combined Regression
lm_combo <- lm(data = df_model_data,
   formula = log10(RR_seq) ~ log10(RR_trips+zero_correction) + 
     log10(RR_move) + 
     log10(RR_air) : flight_length)
lm_combo %>% model_summary()
lm_combo_plot <- plot_pred_RR(lm_combo,y_obs,"Multivariate Model")

#Spline combined regression
slm_combo <- lm(data = df_model_data,
               formula = log10(RR_seq) ~ log10(RR_trips+zero_correction) + 
                 log10(RR_move) + 
                 bs(log10(RR_air), df = 6) * bs(min_cbsa_dist,df=6))
slm_combo %>% model_summary(calc_bootstrap = FALSE)
slm_combo_plot <- plot_pred_RR(slm_combo,y_obs,"Multivariate Model - Spline")

OE_plots <- patchwork::wrap_plots(lm_trips_plot, slm_trips_plot,lm_move_plot, slm_move_plot,
                        lm_air_plot, slm_air_plot,lm_air_int_plot, slm_air_int_plot,
                        lm_euclid_plot, slm_euclid_plot,lm_cbsa_plot, slm_cbsa_plot,
                        patchwork::plot_spacer(),lm_combo_plot, slm_combo_plot,patchwork::plot_spacer(),
                      ncol = 4,nrow=4)
ggsave(paste0("figs/",scenario,"/dist/regression_fit.jpg"),
       plot = OE_plots,
       width = 22,
       height = 20,
       units = "in",
       dpi = 300)
       
