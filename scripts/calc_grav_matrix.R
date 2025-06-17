library(data.table)
library(tidyverse)

df_dist <- fread("data/nb_dist_states.tsv") 
df_pop <- fread("data/state_pop.tsv")

df_grav <- left_join(df_dist,df_pop,by=join_by(state_x==state)) %>%
  rename(pop_x = pop) %>%
  left_join(df_pop,by=join_by(state_y==state)) %>%
  rename(pop_y = pop) %>%
  mutate(grav_wt = ifelse(euclid_dist == 0, NA,
                          as.numeric(pop_x) * as.numeric(pop_y)/(euclid_dist)^2)) %>%
  mutate(log_grav = log10(grav_wt)) %>%
  select(state_x,state_y,euclid_dist,grav_wt,log_grav)

write_tsv(df_grav,"data/state_grav.tsv")
