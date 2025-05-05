library(tidyverse)
library(data.table)
library(dbplyr)
library(duckdb)
library(ggplot2)
library(tictoc)
library(scales)

select <- dplyr::select #Default select to dplyr
SCENARIO <-"USA"

fn_db <- paste0("db_files/db_",SCENARIO,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db,read_only = TRUE)

df_meta <- tbl(con,"metadata")
df_pairs <- tbl(con,"pairs_time")

HYPER_FILTER <- 5E2

strain_counts <- df_pairs %>% 
  group_by(strain_1) %>% 
  summarize(
    n = n(),
    Nextstrain_clade = first(Nextstrain_clade),
  ) %>% 
  filter(n < HYPER_FILTER)

df_pairs_filtered <- df_pairs %>%
  semi_join(strain_counts, by = "strain_1")

state_dict <- df_meta %>% select(strain,division)

#Bind the exposures
df_pairs_filtered <- df_pairs_filtered %>% 
  inner_join(state_dict,join_by(strain_1==strain)) %>%
  rename(x = division) %>% 
  inner_join(state_dict,join_by(strain_2==strain)) %>%
  rename(y = division) %>% 
  filter(!is.na(x) & !is.na(y)) %>%
  filter(x != "NA" & y != "NA")

pair_counts <- df_pairs_filtered %>%
  group_by(x, y) %>%
  summarise(pair_count = n(), .groups = "drop") %>%
  arrange(x) %>%
  collect()

ggplot(pair_counts,aes(x=pair_count)) +
  geom_histogram(fill="royalblue",color="black") +
  scale_x_log10("Number of Identical Pairs Between States",
                labels=scales::label_comma(),
                limits = c(1,1e8)) +
  scale_y_continuous("Number of State Pairs") + 
  labs(title=paste0("Hyper Threshold = ",format(HYPER_FILTER,big.mark = ",",scientific = FALSE))) +
  theme_bw()

ggplot(pair_counts, aes(x = x, y = y, fill = pair_count)) +
  geom_tile() +
  scale_fill_gradient(
    "# of pairs",
    trans = 'log10',
    limits = c(1, 1e8),
    labels = comma_format(),
    low = "darkred",
    high = "lightyellow",
    na.value = "royalblue"
  ) +
  labs(x = NULL, y = NULL,
       title=paste0("Hyper Threshold = ",format(HYPER_FILTER,big.mark = ",",scientific = FALSE))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 5)
  )

