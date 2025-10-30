library(tidyverse)
library(duckdb)
library(dbplyr)
library(zoo)
library(patchwork)
library(lemon)
library(argparse)

source("scripts/cam_map.R")
source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', default = "CAM_1000",
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)
meta_tbl <- tbl(con,'metadata') |> collect()
DBI::dbDisconnect(con)

NATL_POP <- data.frame(
  country = c("Canada","USA","Mexico"),
  natl_pop = c(36991981,331846456,126014024)
)
division_pop <- read_csv("data/division_pop.csv")


# Load shapefile and RR file
cam_map <- prep_cam_map("data/shp-files/cam-shp.gpkg")
df_regions <- read_csv("data/regions.csv") %>%
  select(state,bea_reg)
bea_map <- plot_cam_choropleth(cam_map,
                    data = df_regions,
                    state_col = state,
                    value_col = bea_reg,
                    scale_fun = region_fill_scale,
                    title = "State/Provincial Regions")
ggsave(plot=bea_map,
       filename=paste0("figs/",scenario,"/bea_region_map.jpg"),
       height=7,width=7,dpi=192)

seq_by_div <- meta_tbl |>
  group_by(division) |>
  summarize(num_seq = n()) |>
  select(division,num_seq) |>
  left_join(division_pop,by = join_by(division == division)) |>
  mutate(seq_effort = num_seq / population * 1E5)

seq_map <- plot_cam_choropleth(cam_map,
                    data = seq_by_div,
                    state_col = division,
                    value_col = num_seq,
                    fill_mapper = log10,
                    scale_fun = log_seq_scale,
                    title = "Sequence Counts",
                    bottom_legend = TRUE)   

effort_cap <- function(x,UB = 5000){
  pmin(x,UB)
}

seq_eff_map <- plot_cam_choropleth(cam_map,
                    data = seq_by_div,
                    state_col = division,
                    value_col = seq_effort,
                    fill_mapper = effort_cap,
                    scale_fun = effort_scale,
                    title = "Sequencing per Capita",
                    bottom_legend = TRUE) 

#Effort calculations
week_start_template <- seq.Date(from=as.Date("2020-01-01"),to=as.Date("2024-12-31"),by="week")
empty_effort_template <- tibble(week_start = week_start_template)

natl_week_effort <- meta_tbl |>
  mutate(date = ymd(date)) |>
  filter(!is.na(date)) |>  
  mutate(week_start = floor_date(date, "week")) |>
  group_by(week_start, country) |>
  summarize(num_seq = n()) |>
  select(country, week_start, num_seq) |>
  left_join(NATL_POP,join_by(country == country)) |>
  arrange(week_start) |>
  group_by(country) |>
  mutate(sma_seq = rollmean(num_seq,k=3,fill=NA,align="center")) |>
  mutate(seq_effort = sma_seq/natl_pop * 1E5) |>
  ungroup()
  
division_quarter_effort <- meta_tbl |>
  mutate(date = ymd(date)) |>
  filter(!is.na(date)) |>  
  mutate(quarter = floor_date(date, "quarter")) |>
  group_by(quarter, division) |>
  summarize(num_seq = n(),
            country = unique(country)) |>
  select(division, country, quarter, num_seq) |>
  left_join(division_pop,join_by(division == division)) |>
  arrange(quarter) |>
  group_by(division) |>
  mutate(seq_effort = num_seq/population * 1E5) |>
  ungroup()

plot_natl_effort <- ggplot(natl_week_effort, 
       aes(x = week_start, y = seq_effort, color = country)) +
  geom_line(linewidth = 1.0) +
  scale_x_date(date_breaks = "4 months",name = "Date") +
  scale_y_continuous(name = "Weekly Sequences (per 100,000)") +
  theme_bw()  +
  ggtitle("Country Sequences Per Capita") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  country_color_scale()

plot_effort_box <- ggplot(division_quarter_effort %>% filter(country != "Mexico"), 
       aes(x = quarter, y = seq_effort,fill=country,group = interaction(country, quarter))) +
  geom_boxplot() + 
  theme_bw() +
  ggtitle("State/Provincial Sequences Per Capita") +
  scale_x_date(name = "Date",
               date_breaks = "3 months", 
               limits=c(as.Date("2020-01-01"),
                        as.Date("2024-11-30")),
               date_labels = "%b\n%Y") +
  scale_y_continuous(name = "Quarterly Sequences (per 100,000)",
                     limits = c(0,1000)) +
  country_fill_scale()

age_sex_seq <- meta_tbl |>
  group_by(sex,age_adj) |>
  summarize(num_seq = n()) |>
  select(age_adj,sex,num_seq)

pyramid_limit <- age_sex_seq %>%
  filter(!is.na(sex)) %>%
  filter(!is.na(age_adj)) %>%
  pull(num_seq) %>% max()

plot_age_sex <- ggplot(
  data = age_sex_seq,
  aes(x=age_adj,
      fill = sex,
      y = ifelse(sex == "Male",num_seq,-num_seq))) +
  geom_bar(stat="identity") +
  scale_x_continuous(
    limits = c(0,90),
    name = "Age"
  ) +
  scale_y_continuous(
    labels = abs, 
    limits = pyramid_limit * c(-1,1),
    name = "Number of Sequences"
  ) +
  sex_fill_scale() +
  coord_flip() +
  theme_bw()  + 
  theme(legend.position = "bottom")

sex_seq <- age_sex_seq %>%
  filter(!is.na(sex)) %>%
  group_by(sex) %>%
  summarize(num_seq = sum(num_seq)) %>%
  mutate(perc_sex = round(100 *num_seq/sum(num_seq),1))
  
plot_age_region <- ggplot(
  data = meta_tbl %>%
    mutate(bea_reg = fct_reorder(bea_reg, age_adj, .fun = median)),
  aes(x = age_adj, y = bea_reg, fill = bea_reg)) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width = 0.3, show.legend = FALSE, outlier.shape = NA) + 
  theme_bw() +
  labs(x="Age",y="Region") +
  scale_x_continuous(limits=c(0,100)) +
  region_fill_scale()

reg_clade_prop <- meta_tbl %>%
  group_by(bea_reg) %>%
  count(clade_nextstrain) %>%
  mutate(perc_clade = n/sum(n)) %>%
  select(clade_nextstrain,bea_reg,perc_clade)

ggplot(reg_clade_prop, aes(fill=clade_nextstrain, y=perc_clade, x=bea_reg)) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal() + coord_flip() +
  labs(y = "Proportion", x = "Region",fill="Nextstrain Clades") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 5))


clade_week_prop_reg <- meta_tbl |>
  mutate(date = ymd(date)) |>
  filter(!is.na(date)) |>  
  mutate(week_start = floor_date(date, "week")) |>
  group_by(week_start, bea_reg, clade_nextstrain) |>
  summarize(num_seq = n(), .groups = "drop") |>
  group_by(bea_reg, clade_nextstrain) |>
  arrange(week_start) |>
  mutate(sma_seq = rollmean(num_seq, k=5, fill=0, align="center")) |>
  group_by(bea_reg, week_start) |>  # Group by week for proportions
  mutate(prop_clade = sma_seq/sum(sma_seq, na.rm = TRUE)) |>
  ungroup()


plot_clade_reg_time <- ggplot(clade_week_prop_reg, aes(x = week_start, y = prop_clade, fill=clade_nextstrain)) +
  geom_area(color="lightgray",
            linewidth = 0.7,
            alpha=0.7,
            position = "identity") +
  labs(
    title = "Proportion of SARS-CoV-2 Sequence Counts by Clade",
    x = "Weekly Start Date",
    y = "Proportion of Sequences",
    fill = "Nextstrain Clade"
  ) +
  scale_x_date(date_breaks = "4 months", 
               date_labels = "%b-%y",
               limits = c(as.Date("2020-01-01"),
                          as.Date("2024-12-31"))) +
  theme_minimal() +  
  theme(
    legend.title = element_text(size = rel(0.7)),  
    legend.text = element_text(size = rel(0.6))   
  ) +
  guides(fill = guide_legend(nrow = 5)) +
  facet_wrap(vars(bea_reg), ncol = 2) 


reposition_legend(plot_clade_reg_time,
                  'center',
                  panel="panel-2-6")

#Stitch the main figures together
combined_fig <- seq_map + seq_eff_map + 
  plot_age_sex + (plot_natl_effort / plot_effort_box) + 
  plot_layout(ncol=2,
            byrow = TRUE) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0.02, 0.98)  
  )
ggsave(combined_fig,
       filename=paste0("figs/",scenario,"/figure_1.jpg"),
       units = "in",
       dpi = 192,
       height = 20,
       width = 18)
