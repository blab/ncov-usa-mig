library(tidyverse, quietly = T)
library(dbplyr, quietly = T)
library(duckdb, quietly = T)
library(ggplot2, quietly = T)
library(argparse, quietly = T)
library(lubridate, quietly = T)
library(zoo, quietly = T)
library(data.table,quietly = T)
library(usmap, quietly = T)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

US_POP <- 331449281 #Total US Pop for Census 2020
STATE_POP <- fread("data/state_pop.tsv") |>
  tibble() |>
  filter(!(state %in% c("Puerto Rico","Northern Mariana","Virgin Islands","Guam","American Samoa")))
fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)
meta_tbl <- tbl(con,'metadata') |> collect()

week_start_template <- seq.Date(from=as.Date("2020-01-01"),to=as.Date("2024-12-31"),by="week")
empty_effort_template <- tibble(week_start = week_start_template) #Fills empty holes
#Weeks are slightly off by year, need to align over the entire time period
find_week_start <- function(x){ 
  max(week_start_template[x > week_start_template])
}
find_week_start <- Vectorize(find_week_start)


natl_week_effort <- meta_tbl |>
  mutate(date = ymd(date)) |>
  filter(!is.na(date)) |>
  mutate(week_start = find_week_start(date) |> as.Date()) |>
  group_by(week_start) |>
  summarize(num_seq = n()) |>
  ungroup () |>
  full_join(empty_effort_template, by = join_by(week_start)) |>
  mutate(num_seq = ifelse(is.na(num_seq),0,num_seq)) |>
  mutate(seq_effort = num_seq/US_POP * 1E5) |>
  arrange(week_start) |>
  mutate(
    sma_effort = rollmean(seq_effort, k = 3, fill = NA, align = "center")
  )

natl_week_effort_state <- tibble(
  num_seq = double(),
  week_start = date(),
  seq_effort = double(),
  sma_effort = double(),
  division = character()
)


for(s in STATE_POP$state){
  pop <- filter(STATE_POP,STATE_POP$state == s) |> pull(pop)
  natl_week_effort_state <- meta_tbl |>
    filter(division == s) |>
    mutate(date = ymd(date)) |>
    filter(!is.na(date)) |>
    mutate(week_start = find_week_start(date) |> as.Date()) |>
    group_by(week_start) |>
    summarize(num_seq = n()) |>
    ungroup () |>
    full_join(empty_effort_template, by = join_by(week_start)) |>
    mutate(num_seq = ifelse(is.na(num_seq),0,num_seq)) |>
    mutate(seq_effort = num_seq/pop * 1E5) |>
    arrange(week_start) |>
    mutate(sma_effort = rollmean(seq_effort, k = 3, fill = NA, align = "center")) |>
    mutate(division = s) |>
    rbind(natl_week_effort_state)
}

clade_list <- unique(meta_tbl$Nextstrain_clade)
natl_week_effort_clade <- tibble(
  num_seq = double(),
  week_start = date(),
  sma_seq = double(),
  clade = character()
)

for(c in clade_list){
  natl_week_effort_clade <- meta_tbl |>
    filter(Nextstrain_clade == c) |>
    mutate(date = ymd(date)) |>
    filter(!is.na(date)) |>
    mutate(week_start = find_week_start(date) |> as.Date()) |>
    group_by(week_start) |>
    summarize(num_seq = n()) |>
    ungroup () |>
    full_join(empty_effort_template, by = join_by(week_start)) |>
    mutate(num_seq = ifelse(is.na(num_seq),0,num_seq)) |>
    arrange(week_start) |>
    mutate(sma_seq = rollmean(num_seq, k = 3, fill = NA, align = "center")) |>
    mutate(clade = c) |>
    rbind(natl_week_effort_clade)
}

natl_week_effort_clade <- natl_week_effort_clade |>
  group_by(week_start) |>
  mutate(prop_strain = num_seq/sum(num_seq)) |>
  ungroup()


state_effort <- tibble(
  num_seq = double(),
  seq_effort = double(),
  state = character()
)
for(s in STATE_POP$state){
  pop <- filter(STATE_POP,STATE_POP$state == s) |> pull(pop)
  state_effort <- meta_tbl |>
    filter(division == s) |>
    summarize(num_seq = n()) |>
    mutate(seq_effort = num_seq/pop * 1E5)  |>
    mutate(state = s) |>
    rbind(state_effort)
}

state_effort <- state_effort %>%
  mutate(rel_effort = seq_effort/(sum(num_seq)/US_POP * 1E5))

state_effort_year <- tibble(
  num_seq = double(),
  year = integer(),
  seq_effort = double(),
  state = character()
)
for(s in STATE_POP$state){
  pop <- filter(STATE_POP,STATE_POP$state == s) |> pull(pop)
  state_effort_year <- meta_tbl |>
    filter(division == s) |>
    mutate(date = ymd(date)) |>
    filter(!is.na(date)) |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarize(num_seq = n()) |>
    mutate(seq_effort = num_seq/pop * 1E5)  |>
    ungroup() |>
    mutate(state = s) |>
    rbind(state_effort_year)
}

state_effort_year <- state_effort_year |>
  group_by(year) |> 
  mutate(rel_effort = seq_effort/(sum(num_seq)/US_POP * 1E5))


#Use these for labels
natl_effort_allyears <- sum(state_effort$num_seq)/US_POP * 1E5
natl_effort_2020 <- sum(state_effort_year |> filter(year == 2020) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar
natl_effort_2021 <- sum(state_effort_year |> filter(year == 2021) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar
natl_effort_2022 <- sum(state_effort_year |> filter(year == 2022) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar
natl_effort_2023 <- sum(state_effort_year |> filter(year == 2023) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar
natl_effort_2024 <- sum(state_effort_year |> filter(year == 2024) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar

plot_seq_natl <- ggplot(natl_week_effort, aes(x = week_start, y = seq_effort)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line(aes(x = week_start, y =sma_effort),
            linetype = "dashed",
            linewidth = 0.8) +
  labs(
    title = "Weekly SARS-CoV-2 Sequencing Effort",
    x = "Weekly Start Date",
    y = "Sequences Per 100,000 Persons"
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%y") + 
  theme_bw()
ggsave("figs/effort/natl_seq_plot.jpg",plot_seq_natl,dpi=300,
       width = 10,
       height = 6,
       create.dir=TRUE)

plot_seq_state <- ggplot(natl_week_effort, aes(x = week_start, y = seq_effort)) +
  #geom_point(data=natl_week_effort_state,
  #         aes(week_start,seq_effort,color=division),
  #         alpha = 0.2,size = 0.5) +
  geom_line(data=natl_week_effort_state,
            aes(week_start,sma_effort,color=division),
            alpha = 0.4,
            linetype = "dashed",) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_line(aes(x = week_start, y =sma_effort),
            linetype = "dashed",
            linewidth = 0.8) +
  labs(
    title = "Weekly SARS-CoV-2 Sequencing Effort",
    x = "Weekly Start Date",
    y = "Sequences Per 100,000 Persons"
  ) +
  theme_bw() +
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%y") + 
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=5))
ggsave("figs/effort/state_seq_plot.jpg",plot_seq_state,dpi=300,
       width = 15,
       height = 9)

plot_seq_clade <- ggplot(natl_week_effort_clade |> filter(num_seq > 0), aes(x = week_start, y = sma_seq, color = clade)) +
  geom_point(data=natl_week_effort,aes(x=week_start,y=num_seq), alpha = 0.6, size = 1.5, color="black") +
  geom_area(data=natl_week_effort,aes(x=week_start,y=num_seq),
            linetype = "dashed",
            linewidth = 0.8,
            alpha = 0.1,
            color="black") +
    geom_point(aes(y=num_seq),alpha = 0.6, size = 1.5) +
  geom_line(linetype = "dashed",
            linewidth = 0.8,
            alpha = 0.4) +
  labs(
    title = "Weekly SARS-CoV-2 Sequence Counts by Clade",
    x = "Weekly Start Date",
    y = "Number of Sequences"
  ) +
  theme_bw() +
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%y") + 
  theme(legend.position = "bottom") +
  guides(color=guide_legend(nrow=5))
ggsave("figs/effort/clade_seq_plot.jpg",plot_seq_clade,dpi=300,
       width = 15,
       height = 9)

plot_seq_clade_area <- ggplot(natl_week_effort_clade, aes(x = week_start, y = sma_seq, fill= clade)) +
  geom_area(color="lightgray",
            linewidth = 0.4) +
  labs(
    title = "Weekly SARS-CoV-2 Sequence Counts by Clade",
    x = "Weekly Start Date",
    y = "SMA of Number of Sequences"
  ) +
  theme_bw() +
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%y") + 
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=5))
ggsave("figs/effort/clade_area_plot.jpg",plot_seq_clade_area,dpi=300,
       width = 15,
       height = 9)

plot_seq_clade_prop <- ggplot(natl_week_effort_clade, aes(x = week_start, y = prop_strain, fill= clade)) +
  geom_area(color="lightgray",
            linewidth = 0.4) +
  labs(
    title = "Proportion of SARS-CoV-2 Sequence Counts by Clade",
    x = "Weekly Start Date",
    y = "Proportion of Sequences"
  ) +
  theme_bw() +
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%y") + 
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(nrow=5))
ggsave("figs/effort/clade_prop_plot.jpg",plot_seq_clade_prop,dpi=300,
       width = 15,
       height = 9)

map_state_effort <- plot_usmap(data = state_effort, region="states", values = "seq_effort") +
  labs(
    title = "Sequencing Effort by State - All Years",
    subtitle = paste0("National Effort - ",formatC(natl_effort_allyears,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradientn(limits=c(0,6000),
                        name="Sequences per 100,000 residents",
                       colors = c("dodgerblue2","seagreen","yellow")
                        ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_all.jpg",map_state_effort,dpi=300,width = 10.5, height = 5, units = "in")

map_state_effort_rel <- plot_usmap(data = state_effort, region="states", values = "rel_effort") +
  labs(
    title = "Relative Effort by State - All Years",
    subtitle = paste0("National Effort - ",formatC(natl_effort_allyears,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradient2(limits = c(0.2,4),
                       breaks = c(0.2,0.5,1,2,4),
                       name="Relative Effort",
                       low = "royalblue4",
                       high = "goldenrod2",
                       transform = "log",
                       midpoint = 1
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_all_rel.jpg",map_state_effort_rel,dpi=300,width = 10.5, height = 5, units = "in")


map_state_effort_2020 <- plot_usmap(data = state_effort_year |> filter(year == 2020), 
                                    region="states", values = "seq_effort") +
  labs(
    title = "Sequencing Effort by State - 2020 ",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2020,big.mark = ",",format="d")," sequences per 100,000 residents")
    ) +
  scale_fill_gradientn(limits=c(0,2600),
                       name="Sequences per 100,000 residents",
                       colors = c("dodgerblue2","seagreen","yellow")
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_2020.jpg",map_state_effort_2020,dpi=300,width = 10.5, height = 5, units = "in")

map_state_effort_rel_2020 <- plot_usmap(data = state_effort_year |> filter(year == 2020), region="states", values = "rel_effort") +
  labs(
    title = "Relative Effort by State - 2020",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2020,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradient2(breaks = c(0.05,0.1,0.2,0.5,1,2,5,10,20),
                       name="Relative Effort",
                       low = "royalblue4",
                       high = "goldenrod2",
                       transform = "log",
                       midpoint = 1
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_all_rel_2020.jpg",map_state_effort_rel_2020,dpi=300,width = 10.5, height = 5, units = "in")

natl_effort_2021 <- sum(state_effort_year |> filter(year == 2021) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar
map_state_effort_2021 <- plot_usmap(data = state_effort_year |> filter(year == 2021), 
                                    region="states", values = "seq_effort") +
  labs(
    title = "Sequencing Effort by State - 2021",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2021,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradientn(limits=c(0,2600),
                       name="Sequences per 100,000 residents",
                       colors = c("dodgerblue2","seagreen","yellow")
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_2021.jpg",map_state_effort_2021,dpi=300,width = 10.5, height = 5, units = "in")

map_state_effort_rel_2021 <- plot_usmap(data = state_effort_year |> filter(year == 2021), region="states", values = "rel_effort") +
  labs(
    title = "Relative Effort by State - 2021",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2021,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradient2(breaks = c(0.05,0.1,0.2,0.5,1,2,5,10,20),
                       name="Relative Effort",
                       low = "royalblue4",
                       high = "goldenrod2",
                       transform = "log",
                       midpoint = 1
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_all_rel_2021.jpg",map_state_effort_rel_2021,dpi=300,width = 10.5, height = 5, units = "in")

natl_effort_2022 <- sum(state_effort_year |> filter(year == 2022) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar
map_state_effort_2022 <- plot_usmap(data = state_effort_year |> filter(year == 2022), 
                                    region="states", values = "seq_effort") +
  labs(
    title = "Sequencing Effort by State - 2022",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2022,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradientn(limits=c(0,2600),
                       name="Sequences per 100,000 residents",
                       colors = c("dodgerblue2","seagreen","yellow")
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_2022.jpg",map_state_effort_2022,dpi=300,width = 10.5, height = 5, units = "in")

map_state_effort_rel_2022 <- plot_usmap(data = state_effort_year |> filter(year == 2022), region="states", values = "rel_effort") +
  labs(
    title = "Relative Effort by State - 2022",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2022,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradient2(breaks = c(0.05,0.1,0.2,0.5,1,2,5,10,20),
                       name="Relative Effort",
                       low = "royalblue4",
                       high = "goldenrod2",
                       transform = "log",
                       midpoint = 1
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_all_rel_2022.jpg",map_state_effort_rel_2022,dpi=300,width = 10.5, height = 5, units = "in")

natl_effort_2023 <- sum(state_effort_year |> filter(year == 2023) |> pull(num_seq))/US_POP * 1E5 #Excuse my grammar
map_state_effort_2023 <- plot_usmap(data = state_effort_year |> filter(year == 2023), 
                                    region="states", values = "seq_effort") +
  labs(
    title = "Sequencing Effort by State - 2023",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2023,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradientn(limits=c(0,2600),
                       name="Sequences per 100,000 residents",
                       colors = c("dodgerblue2","seagreen","yellow")
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_2023.jpg",map_state_effort_2023,dpi=300,width = 10.5, height = 5, units = "in")

map_state_effort_rel_2023 <- plot_usmap(data = state_effort_year |> filter(year == 2023), region="states", values = "rel_effort") +
  labs(
    title = "Relative Effort by State - 2023",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2023,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradient2(breaks = c(0.05,0.1,0.2,0.5,1,2,5,10,20),
                       name="Relative Effort",
                       low = "royalblue4",
                       high = "goldenrod2",
                       transform = "log",
                       midpoint = 1
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_all_rel_2023.jpg",map_state_effort_rel_2023,dpi=300,width = 10.5, height = 5, units = "in")

map_state_effort_2024 <- plot_usmap(data = state_effort_year |> filter(year == 2024), 
                                    region="states", values = "seq_effort") +
  labs(
    title = "Sequencing Effort by State - 2024",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2024,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradientn(limits=c(0,200),
                       name="Sequences per 100,000 residents",
                       colors = c("dodgerblue2","seagreen","yellow")
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_2024.jpg",map_state_effort_2024,dpi=300,width = 10.5, height = 5, units = "in")

map_state_effort_rel_2024 <- plot_usmap(data = state_effort_year |> filter(year == 2024), region="states", values = "rel_effort") +
  labs(
    title = "Relative Effort by State - 2024",
    subtitle = paste0("National Effort - ",formatC(natl_effort_2024,big.mark = ",",format="d")," sequences per 100,000 residents")
  ) +
  scale_fill_gradient2(breaks = c(0.05,0.1,0.2,0.5,1,2,5,10,20),
                       name="Relative Effort",
                       low = "royalblue4",
                       high = "goldenrod2",
                       transform = "log",
                       midpoint = 1
  ) +
  theme(legend.position=c(0.9,0.05))
ggsave("figs/effort/map_state_all_rel_2024.jpg",map_state_effort_rel_2024,dpi=300,width = 10.5, height = 5, units = "in")

ggplot(data=state_effort_year |> arrange(rel_effort),aes(x=year,y=rel_effort,color=state)) +
  geom_point(alpha = 0.1) +
  geom_line(alpha = 0.2) +
  scale_y_continuous(transform="log",
                     breaks = c(0.03, 0.1, 0.5, 1, 2, 10, 30)) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Relative Effort"
  )

write.csv(natl_week_effort,"results/USA/natl_effort.csv")
