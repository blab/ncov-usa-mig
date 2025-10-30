library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(argparse)
source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
SCENARIO <- args$scenario
fig_path <- paste0("figs/",SCENARIO,"/time/")

state_rr_snap <- read_tsv(paste0("results/",
                                 SCENARIO,
                                 "/time_state/df_state_rr_snap.tsv"))
state_rr_series <- read_tsv(paste0("results/",
                                   SCENARIO,
                                   "/time_state/df_state_rr_series.tsv"))
state_pair_types <- state_rr_series %>%
  select(x, y, pair_type) %>%
  distinct()
state_rr_snap <- state_rr_snap %>%
  left_join(state_pair_types, by = c("x", "y"))

df_regions <- read_csv("data/us_states_regions.csv") %>%
  select(state, bea_reg)


PAIR_LIST <- data.frame(
  x = c("Mexico",
        "Mexico",
        "British Columbia",
        "Mexico",
        "Ontario",
        "Ontario",
        "Ontario",
        "Ontario",
        "Quebec",
        "California",
        "California",
        "California",
        "California",
        "Illinois"),
  y = c("Texas",
        "California",
        "Washington",
        "Arizona",
        "Michigan",
        "New York",
        "Quebec",
        "British Columbia",
        "New York",
        "New York",
        "Arizona",
        "Texas",
        "Washington",
        "Michigan")
)

#Make quarter labels
format_quarters <- function(date){
  paste0(lubridate::year(date)," Q",lubridate::quarter(date))
}

sub_RR_snap <- state_rr_snap %>%
  inner_join(PAIR_LIST, by=c("x","y")) %>%
  mutate(pair_name = paste0(x,"-",y)) %>%
  arrange(pair_name) %>%
  mutate(is_intl = pair_type == "International") %>%
  # Create separate factor levels for each facet
  group_by(is_intl) %>%
  mutate(pair_name = factor(pair_name, levels = unique(pair_name[order(pair_name)]))) %>%
  mutate(quarter_label = format_quarters(date)) %>%
  # Add fold change calculation
  group_by(pair_name) %>%
  arrange(date) %>%
  mutate(nRR_fold = nRR / first(nRR)) %>%
  ungroup()

# Create common scale for both plots
common_fill_scale <- scale_fill_distiller(
  palette = "Spectral",
  limits = c(-1, 1),
  oob = scales::squish,
  name = "Fold Change",
  labels = function(x) format(10^x, digits = 1)
)

# Create separate plots for each facet
p1 <- ggplot(sub_RR_snap %>% filter(!is_intl),
       aes(x=date, y=pair_name, fill = log10(nRR_fold))) +
  geom_tile(color="#555") +
  scale_x_date(
    date_breaks = "3 months",
    labels = function(x) format_quarters(x),
    expand = expansion(0.01)
  ) +
  labs(x = "Date", y = "Domestic Pairs") +
  common_fill_scale +  # Use common scale
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.justification = "center",
    panel.grid.major = element_blank()
  )

p2 <- ggplot(sub_RR_snap %>% filter(is_intl),
       aes(x=date, y=pair_name, fill=log10(nRR_fold))) +
  geom_tile(color="#555") +
  scale_x_date(
    date_breaks = "3 months",
    labels = function(x) format_quarters(x),
    expand = expansion(0.01)
  ) +
  labs(x = "Date", y = "International Pairs") +
  common_fill_scale +  # Use common scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size=8),
    legend.position = "none",
    panel.grid.major = element_blank()
  ) 
# Combine plots using patchwork
p1 / p2 +
  plot_annotation(
    title = "Normalized RR Over Time by State/Province Pair",
    theme = theme_minimal())

ggsave(paste0(fig_path,"state_pair_nRR_heatmap.png"),
       width=8,
       height=5,
       dpi=192)

unique_series <- state_rr_series %>% 
  mutate(duplicate = x >= y) %>% 
  filter (duplicate != TRUE) %>%
  mutate(pair_name = paste0(x,"-",y)) %>%
  arrange(pair_name) %>%
  mutate(pair_type = factor(pair_type,
                            levels = c("Same State/Province",
                                       "Different States, Same Region",
                                       "Different Regions, Same Country",
                                       "International")))

ggplot(unique_series,
       aes(x = date, y = nRR, group = date, colour = pair_type)) +
  geom_boxplot(alpha = 0.2) +
  geom_smooth(aes(group = 1)) +
  coord_cartesian(ylim=c(0,0.5)) +
  facet_wrap(~ pair_type, ncol = 1) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(paste0(fig_path,"rr_series.png"),
       width=8,
       height=8,
       dpi=192)

ts_nRR_df <- unique_series %>%
  mutate(
    date_num  = as.numeric(date),
    pair_type = factor(pair_type),
    pair_name = factor(pair_name)
  ) %>%
  filter(!is.na(nRR), !is.na(date_num), !is.na(pair_type), !is.na(pair_name))


#ts_nRR_model <- mgcv::gam(
#  nRR ~ s(date_num,bs="cs",k=10) + pair_type + s(pair_name,bs="re"),
#  data=ts_nRR_df
#)

#summary(ts_nRR_model)
#mgcv::gam.check(ts_nRR_model)
#AIC(ts_nRR_model)

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

state_list <- list("California",
                   "Washington",
                   "New York",
                   "New Jersey",
                   "Illinois",
                   "Wisconsin",
                   "Mexico",
                   "Ontario")

state_nrr_timeline <- function(state_name){
    nrr_series_sub <- state_rr_series %>%
        filter(x == state_name & y != state_name) %>%
        rename(paired_state = y) %>%
        select(date, paired_state, nRR) %>%
        left_join(df_regions, by = c("paired_state" = "state"))
        
    plot <- ggplot(nrr_series_sub,
                   aes(x = date, y = nRR, 
                       group = paired_state, 
                       color = bea_reg)) +
        geom_line(alpha = 0.5, linewidth = 0.3) +
        geom_smooth(aes(group=bea_reg), size=1) +
        scale_x_date(name = "Date",date_breaks = "6 months") +
        scale_y_continuous(name = "nRR") +
        ggtitle(paste0("nRR Time Series - ",state_name)) +
        region_color_scale() +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    
    fn_plot <- paste0(fig_path,
                      "/individual_states/nRR_timeline_",
                      state_name,
                      ".png")
    ggsave(plot = plot,
        filename=fn_plot,
        width=8,
        height=5,
        dpi=192,
        create.dir = TRUE)
    return(plot)
}

lapply(state_list,state_nrr_timeline)
