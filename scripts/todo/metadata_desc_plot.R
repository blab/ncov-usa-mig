# Filename: metadata_desc_plots.R
# Author(s): Amin Bemanian
# Date: 07/08/25
# Description: Make descriptive plots for the sequence
# metadata. 
# Arguments: scenario
# TODO: Standardize the color palettes 

library(tidyverse)
library(dbplyr)
library(duckdb)
library(ggpubr)

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

#Debugging shortcuts
scenario <- "CAM_1000"
#ci_flag <- FALSE

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)
df_meta <- tbl(con,"metadata")
df_meta %>%
  count(bea_reg) %>%
  mutate(perc = n / sum(n)) %>%
  collect() %>% arrange(-perc)
DBI::dbDisconnect(con)

#For graphiing
df_meta_full <- df_meta %>% collect()
ggplot(
  data = df_meta_full %>%
    mutate(bea_reg = fct_reorder(bea_reg, age_adj, .fun = median)),
  aes(x = age_adj, y = bea_reg, fill = bea_reg)
) +
  geom_violin(alpha = 0.5, show.legend = FALSE) +
  geom_boxplot(width = 0.3, show.legend = FALSE, outlier.shape = NA) + 
  # stat_compare_means(
  #   method = "wilcox.test",      # or "t.test"
  #   label = "p.signif",          # show stars like **, *** etc.
  #   comparisons = list(c("Eastern Canada", "Western Canada"),
  #                      c("Western Canada", "Mexico"),
  #                      c("Mexico","Southeast"),
  #                      c("Southeast","Mideast"),
  #                      c("Mideast","Plains"),
  #                      c("Plains","Great Lakes"),
  #                      c("Great Lakes", "Southwest"),
  #                      c("Southwest", "New England"),
  #                      c("New England", "Rocky Mountain"),
  #                      c("Rocky Mountain", "Far West")),
  #   hide.ns = TRUE               # optionally hide non-significant
  # ) +
  theme_bw() +
  labs(x="Age",y="Region")

df_clade <- df_meta_full %>%
  mutate(clade_who = ifelse(is.na(clade_who),"Wildtype",clade_who)) %>%
  group_by(bea_reg) %>%
  count(clade_who) %>%
  mutate(perc_clade = n/sum(n))

ggplot(df_clade, aes(fill=clade_who, y=perc_clade, x=bea_reg)) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal() + coord_flip() +
  labs(y = "Proportion", x = "Region",fill="WHO Naming Convention")

df_year <- df_meta_full %>%
  mutate(year = date %>% as.Date() %>% format('%Y')) %>%
  filter(!is.na(year) & year != "2025") %>%
  group_by(bea_reg) %>%
  count(year) %>%
  mutate(perc_year = n/sum(n))

ggplot(df_year, aes(fill=year, y=perc_year, x=bea_reg)) + 
  geom_bar(position="fill", stat="identity") +
  theme_minimal() + coord_flip() +
  labs(y = "Proportion", x= "Region",fill="Year")
