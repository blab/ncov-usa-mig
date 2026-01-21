#File: age_gender_analysis.R
#Author(s): Amin Bemanian
#Date: 12/17/25
#Description: Calculate identical sequence RR across age groups for gender discordant pairs
# This is to help make sure the same age RR patterns we are see aren't due to repeated sampling of the same individuals 
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)
#--ci: TRUE or FALSE for confidence interval calculation

library(tidyverse)
library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(duckdb)
library(dbplyr)
filter <- dplyr::filter

source("scripts/calculate_rr_matrix.R")
source("scripts/bind_pairs_exp.R")
source("scripts/calculate_rr_ci.R")
source("scripts/color_schemes.R")

AGE_UB <- "80y" #Censor older ages group

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', default = "CAM_1000", type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--exclude_duplicates', type = 'logical', default = TRUE, help = "Whether to exclude possible duplicate pairs, default is TRUE")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
exclude_duplicates <- args$exclude_duplicates

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)

gender_pairs <- con %>%
  bind_pairs_exp("sex", exclude_duplicates = exclude_duplicates) %>%
  rename(gender.x = x) %>%
  rename(gender.y = y)

age_pairs <- con %>%
  bind_pairs_exp("age_class", exclude_duplicates = exclude_duplicates)

age_gender_pairs <- age_pairs %>% 
  left_join(gender_pairs)

rr_age_gender_different <- age_gender_pairs %>%
  filter(!is.na(gender.x)) %>%
  filter(!is.na(gender.y)) %>%
  filter(gender.x != gender.y) %>%
  calculate_rr_matrix() %>%
  mutate(same_gender = FALSE)

rr_age_gender_same <- age_gender_pairs %>%
  filter(!is.na(gender.x)) %>%
  filter(!is.na(gender.y)) %>%
  filter(gender.x == gender.y) %>%
  calculate_rr_matrix() %>%
  mutate(same_gender = TRUE)

rr_age_gender <- rbind(
  rr_age_gender_different,
  rr_age_gender_same
) %>% mutate(age_censor = (x > AGE_UB | y > AGE_UB))

NTH_BREAK <- 5
age_levels <- age_pairs %>% 
  pull(x) %>% 
  factor() %>% 
  levels()
age_breaks <- age_levels[seq(1,length(age_levels),by=NTH_BREAK)] 

fig_path <- paste0("figs/",scenario,"/age_gender/")

age_same_plot_by_gender  <- rr_age_gender  %>%
  filter(x==y) %>% #Only look at RR within the same 
  filter(!age_censor) %>%
  ggplot(aes(x = x, colour = same_gender)) +
  geom_point(aes(y = RR)) +
  geom_line(aes(y = RR, group = same_gender)) +
  #geom_linerange(aes(ymin = ci_lb, ymax = ci_ub,group = sameState)) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(name = 'Age Group', breaks = age_breaks) +
  scale_y_continuous(transform ='log',
                     name=expression(RR["identical sequences"]),
                     breaks = c(1E-1, 8E-1, 1, 5, 2, 1E1, 1E2),
                     labels = c(0.1,0.8,1,5,2,10,100),
                     expand = expansion(mult = c(0.18, 0.13)),
                     limits = c(10^(-0.1),10^(1.8))) +
  scale_color_brewer(type="qual", name = "Same Gender",palette = "Dark2") + 
  theme_classic() +
  ggtitle(paste0("Same Age RR")) +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1.),
        strip.text = element_text(size = 12),
        strip.background = element_rect(colour = 'white'),
        plot.title=element_text(hjust=0.5),
        legend.position = "bottom"
      )

ggsave(plot = age_same_plot_by_gender, 
  filename = paste0(fig_path,"same_rr_age_gender.png"),
  dpi = 300,
  height = 5
)

rr_age_gender_equal_age <- rr_age_gender %>%
  filter(x == y) %>%
  mutate(age_cat = case_when(
    x <= "22y" ~ "0-22 yo",
    x <= "44" ~ "23-44 yo",
    x <= "64y" ~ "45-64 yo",
    x <= "74y" ~ "65-74 yo",
    .default = "75+ yo"
  )) %>%
  mutate(age_cat = 
    factor(age_cat,levels = c("0-22 yo",
      "23-44 yo",
      "45-64 yo",
      "65-74 yo",
      "75+ yo"
    )))

same_age_boxplot <- ggplot(rr_age_gender_equal_age,aes(x=age_cat,y=RR,fill=same_gender)) + 
  geom_boxplot() +
  geom_abline(slope = 0, intercept = 0, lty = "dashed") +
  theme_bw() +
  scale_y_log10("Same Age RR",limits=c(0.8,100)) +
  scale_x_discrete("Age Category") +
  scale_fill_brewer(type="qual", name = "Same Gender",palette = "Dark2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "bottom")

ggsave(plot = same_age_boxplot, 
  filename = paste0(fig_path,"same_age_gender_boxplot.png"),
  dpi = 300,
  width = 8,
  units = "in"
)

rr_age_gender_wide <- left_join(
  rr_age_gender_same %>% rename(RR_same = RR) %>% select(x,y,RR_same),
  rr_age_gender_different %>% rename(RR_different = RR) %>% select(x,y,RR_different)
) %>% mutate(age_censor = (x > AGE_UB | y > AGE_UB))

corr_plot_lb <- 0.5
corr_plot_ub <- 30

corr_age_gender <- cor(x=rr_age_gender_wide$RR_same, y=rr_age_gender_wide$RR_different,method = "pearson")

age_gender_corr_plot <- ggplot(rr_age_gender_wide %>% filter(!age_censor), aes(x = RR_same, y = RR_different)) +
  geom_point(alpha = 0.1,color="#36f") +
  geom_abline(slope=1,intercept=0,lty="dashed") +
  scale_y_log10(name="Different Gender RR") +
  scale_x_log10(name="Same Gender RR") +
  labs(title = "Age RR Correlation by Gender Concordance",
    subtitle = sprintf("Pearson's r: %.3f",corr_age_gender)) +
  theme_bw() +
  coord_fixed(ratio = 1,xlim = c(corr_plot_lb,corr_plot_ub), ylim = c(corr_plot_lb,corr_plot_ub))

ggsave(plot = age_gender_corr_plot, 
  filename = paste0(fig_path,"rr_age_gender_corr.png"),
  dpi = 300,
  height = 5
)


#Make heatmaps by gender concordance
UB <- 1.5
LB <- 0.7
fill_bound <- function(x){
  log_x <- log10(x)
  bound_x <- max(min(log_x,log10(UB)),log10(LB))
  return(bound_x)
}
AXIS_SIZE <- 8
make_age_heatmap <- function(data, same_flag, title = "",
                             show_legend = TRUE, show_axis = FALSE, show_rr_labels = FALSE) {
  # Start with age filter
  filtered_data <- data %>%
    filter(!age_censor) %>% 
    filter(same_gender == same_flag)
  
  # Apply fill transformation
  filtered_data <- filtered_data %>%
    rowwise() %>%
    mutate(fill_RR = fill_bound(RR))

  # Create base plot
  p <- ggplot(filtered_data, aes(x=x, y=y, fill=fill_RR)) +
    geom_tile() +
    RR_log_grad(LB = LB,UB = UB) +
    theme_minimal() +
    labs(title = title) +
    theme(plot.title = element_text(hjust=0.5)) +
    coord_equal()

  # Add RR labels if requested
  if (show_rr_labels) {
    p <- p + geom_text_repel(aes(label=round(RR, digits=2)),
                             color = "black",
                             size = RR_SIZE,
                             bg.color = "white",
                             bg.r = 0.1,
                             force = 0)
  }

  # Configure axis display
  if (show_axis) {
    # Get unique age values and filter to every 10 years
    age_breaks <- unique(c(filtered_data$x, filtered_data$y)) %>%
      sort() %>%
      grep("*0y$", ., value = TRUE)
    print(age_breaks)
    p <- p +
      scale_x_discrete(name = "Age Groups", breaks = age_breaks) +
      scale_y_discrete(name = "Age Groups", breaks = age_breaks) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = AXIS_SIZE),
            axis.text.y = element_text(size = AXIS_SIZE))
  } else {
    p <- p + theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  }

  # Optionally hide legend
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

same_gender_age_heatmap <- make_age_heatmap(rr_age_gender,same_flag = TRUE,title = "Same Gender Pairs",show_axis = TRUE)
different_gender_age_heatmap <- make_age_heatmap(rr_age_gender,same_flag = FALSE,title = "Different Gender Pairs",show_axis = TRUE)

ggsave(plot=same_gender_age_heatmap,
  filename=paste0(fig_path,"same_gender_age_heatmap.png"),
  dpi = 300,
  width = 6
)

ggsave(plot=different_gender_age_heatmap,
  filename=paste0(fig_path,"different_gender_age_heatmap.png"),
  dpi = 300,
  width = 6
)

####GENDER RR ANALYSIS
# Function to perform chi-squared test and return formatted string
calc_chisq_label <- function(rr_data) {
  rr_collected <- rr_data %>% collect()

  # Create contingency table from pair counts
  contingency_matrix <- matrix(0, nrow = length(unique(rr_collected$x)),
                                ncol = length(unique(rr_collected$y)))
  rownames(contingency_matrix) <- unique(rr_collected$x)
  colnames(contingency_matrix) <- unique(rr_collected$y)

  for(i in 1:nrow(rr_collected)) {
    row_idx <- which(rownames(contingency_matrix) == rr_collected$x[i])
    col_idx <- which(colnames(contingency_matrix) == rr_collected$y[i])
    contingency_matrix[row_idx, col_idx] <- rr_collected$N_pairs[i]
  }

  chisq_result <- tryCatch({
    chisq.test(contingency_matrix)
  }, error = function(e) {
    return(NULL)
  })

  if(is.null(chisq_result)) {
    return("")
  } else {
    return(sprintf("χ² = %.1f, p %s",
                   chisq_result$statistic,
                   ifelse(chisq_result$p.value < 0.001, "< 0.001",
                          sprintf("= %.3f", chisq_result$p.value))))
  }
}

rr_gender <- con %>%
  bind_pairs_exp("sex",exclude_duplicates = exclude_duplicates) %>%
  calculate_rr_matrix()

gender_pairs_over_age <- gender_pairs %>%
  rename(x = gender.x) %>% rename(y = gender.y)  %>%
  left_join(age_pairs %>% rename(age.x = x) %>% rename(age.y = y)) %>%
  mutate(same_age = (age.x == age.y))

rr_gender_same_age <- gender_pairs_over_age %>%
  filter(same_age) %>%
  calculate_rr_matrix()

rr_gender_diff_age <- gender_pairs_over_age %>%
  filter(!same_age) %>%
  calculate_rr_matrix()

# Add RR matrix for ages off by exactly 1 year (to avoid duplicate issues)
# Use DuckDB SQL functions for string manipulation
# Need to quote column names with dots in them
gender_pairs_over_age_numeric <- gender_pairs_over_age %>%
  mutate(age_x_numeric = sql('CAST(regexp_replace("age.x", \'y.*\', \'\') AS INTEGER)')) %>%
  mutate(age_y_numeric = sql('CAST(regexp_replace("age.y", \'y.*\', \'\') AS INTEGER)')) %>%
  mutate(age_diff = abs(age_x_numeric - age_y_numeric))

rr_gender_age_off_by_1 <- gender_pairs_over_age_numeric %>%
  filter(age_diff == 1) %>%
  select(-age_x_numeric, -age_y_numeric, -age_diff, -age.x, -age.y, -same_age) %>%
  calculate_rr_matrix()

# Calculate chi-squared statistics for each group
chisq_all <- calc_chisq_label(rr_gender)
chisq_same <- calc_chisq_label(rr_gender_same_age)
chisq_diff <- calc_chisq_label(rr_gender_diff_age)
chisq_off1 <- calc_chisq_label(rr_gender_age_off_by_1)

rr_gender_combined <- rbind(
  rr_gender %>% mutate(age_group = paste0("All Pairs\n", chisq_all)),
  rr_gender_same_age %>% mutate(age_group = paste0("Same Age Pairs\n", chisq_same)),
  rr_gender_diff_age %>% mutate(age_group = paste0("Different Age Pairs\n", chisq_diff)),
  rr_gender_age_off_by_1 %>% mutate(age_group = paste0("Age Off by 1 Year\n", chisq_off1))) %>%
  rowwise() %>%
  mutate(fill_RR = fill_bound(RR)) %>%
  ungroup() %>%
  mutate(age_group = factor(age_group, levels = c(
    paste0("All Pairs\n", chisq_all),
    paste0("Different Age Pairs\n", chisq_diff),
    paste0("Age Off by 1 Year\n", chisq_off1),
    paste0("Same Age Pairs\n", chisq_same)
  )))

rr_gender_heatmap <- ggplot(rr_gender_combined, aes(x=x, y=y, fill=fill_RR)) +
  geom_tile() +
  RR_log_grad(LB = LB,UB = UB) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5)) +
  coord_equal() +
  labs(x="",y="",title="Gender RR by Age") +
  ggrepel::geom_text_repel(aes(label=round(RR, digits=2)),
    color = "black",
    size = 4,
    bg.color = "white",
    bg.r = 0.1,
    force = 0) +
  facet_wrap(vars(age_group),nrow=2)

ggsave(plot=rr_gender_heatmap,
  filename = paste0(fig_path,"gender_rr_over_age.png"),
  dpi = 300,
  width = 7.5,
  height = 5,
  units = "in"
)

#### GENDER RR BY LIFE STAGE
# Add life stage classifications in database and calculate RR for each combination
gender_pairs_with_lifestage <- gender_pairs_over_age %>%
  mutate(lifestage.x = sql("CASE
    WHEN CAST(regexp_replace(\"age.x\", 'y.*', '') AS INTEGER) <= 5 THEN 'Early Childhood (0-5)'
    WHEN CAST(regexp_replace(\"age.x\", 'y.*', '') AS INTEGER) <= 11 THEN 'Young Schoolage (6-11)'
    WHEN CAST(regexp_replace(\"age.x\", 'y.*', '') AS INTEGER) <= 17 THEN 'Adolescence (12-17)'
    WHEN CAST(regexp_replace(\"age.x\", 'y.*', '') AS INTEGER) <= 22 THEN 'College (18-22)'
    WHEN CAST(regexp_replace(\"age.x\", 'y.*', '') AS INTEGER) <= 49 THEN 'Early Adulthood (23-49)'
    WHEN CAST(regexp_replace(\"age.x\", 'y.*', '') AS INTEGER) <= 64 THEN 'Middle Adulthood (50-64)'
    WHEN CAST(regexp_replace(\"age.x\", 'y.*', '') AS INTEGER) <= 79 THEN 'Early Seniors (65-79)'
    ELSE 'Late Seniors (80+)' END")) %>%
  mutate(lifestage.y = sql("CASE
    WHEN CAST(regexp_replace(\"age.y\", 'y.*', '') AS INTEGER) <= 5 THEN 'Early Childhood (0-5)'
    WHEN CAST(regexp_replace(\"age.y\", 'y.*', '') AS INTEGER) <= 11 THEN 'Young Schoolage (6-11)'
    WHEN CAST(regexp_replace(\"age.y\", 'y.*', '') AS INTEGER) <= 17 THEN 'Adolescence (12-17)'
    WHEN CAST(regexp_replace(\"age.y\", 'y.*', '') AS INTEGER) <= 22 THEN 'College (18-22)'
    WHEN CAST(regexp_replace(\"age.y\", 'y.*', '') AS INTEGER) <= 49 THEN 'Early Adulthood (23-49)'
    WHEN CAST(regexp_replace(\"age.y\", 'y.*', '') AS INTEGER) <= 64 THEN 'Middle Adulthood (50-64)'
    WHEN CAST(regexp_replace(\"age.y\", 'y.*', '') AS INTEGER) <= 79 THEN 'Early Seniors (65-79)'
    ELSE 'Late Seniors (80+)' END"))

# Define life stage levels for ordering
lifestage_levels <- c(
  "Early Childhood (0-5)",
  "Young Schoolage (6-11)",
  "Adolescence (12-17)",
  "College (18-22)",
  "Early Adulthood (23-49)",
  "Middle Adulthood (50-64)",
  "Early Seniors (65-79)",
  "Late Seniors (80+)"
)

# Calculate RR for each life stage combination
rr_gender_by_lifestage <- NULL
for(ls_x in lifestage_levels) {
  for(ls_y in lifestage_levels) {
    rr_temp <- gender_pairs_with_lifestage %>%
      filter(lifestage.x == ls_x, lifestage.y == ls_y) %>%
      select(strain_1, strain_2, x, y, n_mutations, possible_duplicates) %>%
      calculate_rr_matrix() %>%
      collect() %>%
      mutate(lifestage.x = ls_x, lifestage.y = ls_y)

    rr_gender_by_lifestage <- rbind(rr_gender_by_lifestage, rr_temp)
  }
}

# Apply fill bounds and set factor levels
# Filter to only same life stage pairs (diagonal)
rr_gender_by_lifestage <- rr_gender_by_lifestage %>%
  filter(lifestage.x == lifestage.y) %>%
  rowwise() %>%
  mutate(fill_RR = fill_bound(RR)) %>%
  ungroup() %>%
  mutate(
    lifestage.x = factor(lifestage.x, levels = lifestage_levels)
  )

# Create faceted heatmap by life stage (only within-stage comparisons)
lifestage_gender_heatmap <- ggplot(rr_gender_by_lifestage, aes(x=x, y=y, fill=fill_RR)) +
  geom_tile() +
  RR_log_grad(LB = LB, UB = UB) +
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5),
        strip.text = element_text(size = 9),
        axis.text = element_text(size = 8)) +
  coord_equal() +
  labs(x="", y="", title="Within Life Stage Gender RR") +
  geom_text(aes(label=round(RR, digits=2)),
            color = "black",
            size = 3.5) +
  facet_wrap(vars(lifestage.x), nrow = 2)

ggsave(plot=lifestage_gender_heatmap,
  filename = paste0(fig_path,"gender_rr_by_lifestage.png"),
  dpi = 300,
  width = 12,
  height = 9,
  units = "in"
)

#Disconnect DB safely
DBI::dbDisconnect(con)
