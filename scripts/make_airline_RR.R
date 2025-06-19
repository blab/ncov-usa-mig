library(tidyverse)

EXCL_STATES <- c("PR","TT","VI")

abbrev_to_state <- function(abbrev) {
  lookup <- c(state.abb, "DC")
  names <- c(state.name, "District of Columbia")
  names[match(toupper(abbrev), lookup)]
}

t100_set<-read_csv("data/T100_Dataset/T_T100D_MARKET_ALL_CARRIER.csv") #T100 data, has individual legs (i.e. layovers are reported as trips)
t100_agg  <- t100_set %>%
  filter(!(ORIGIN_STATE_ABR %in% EXCL_STATES)) %>%
  filter(!(DEST_STATE_ABR %in% EXCL_STATES)) %>%
  select(PASSENGERS,DISTANCE,ORIGIN_STATE_ABR,DEST_STATE_ABR,MONTH) %>%
  group_by(ORIGIN_STATE_ABR,DEST_STATE_ABR) %>%
  summarize(pass = sum(PASSENGERS)) %>%
  rename(ori = ORIGIN_STATE_ABR) %>%
  rename(dest = DEST_STATE_ABR) %>%
  ungroup()

uniq_states <- union(t100_agg$ori,t100_agg$dest)
all_grid <- expand_grid(
  ori = uniq_states,
  dest = uniq_states,
  plus_one = 1
)

t100_agg<-right_join(t100_agg,all_grid,
                    by=join_by(ori,dest)) %>%
  mutate(pass = ifelse(is.na(pass),0,pass) + plus_one) %>% #Basically forces every pair to have 1 trip at minimum for RR calculation
  select(!plus_one)
t100_mirror <- t100_agg %>%
  rename(ori_mirror = dest) %>%
  rename(dest_mirror = ori) %>%
  rename(pass_mirror = pass)

df_rr_t100 <- left_join(t100_agg,t100_mirror,
          by = join_by("ori" == "ori_mirror",
                       "dest" == "dest_mirror")) %>%
  mutate(x = abbrev_to_state(ori)) %>%
  mutate(y = abbrev_to_state(dest)) %>%
  mutate(pass_xy = pass + pass_mirror,
         pass_total = sum(pass_xy)) %>%
  group_by(x) %>%
  mutate(pass_x = sum(pass_xy)) %>% ungroup() %>%
  group_by(y) %>%
  mutate(pass_y = sum(pass_xy)) %>% ungroup() %>%
  mutate(RR_air = pass_xy * pass_total / pass_x / pass_y) %>%
  select(x,y,RR_air)
  
write_csv(df_rr_t100,"data/rr_air_t100.csv")

#Sample of 10% of all flight itineraries
t_DB1B_set <- read_csv("data/T100_Dataset/T_DB1B_MARKET.csv") #For DB1B market data, all legs are combined together so its just the full trip OD (no layovers)
t_DB1B_agg <- t_DB1B_set %>%
  filter(!(ORIGIN_STATE_ABR %in% EXCL_STATES)) %>%
  filter(!(DEST_STATE_ABR %in% EXCL_STATES)) %>%
  group_by(ORIGIN_STATE_NM,DEST_STATE_NM) %>%
  summarise(pass = sum(PASSENGERS)) %>%
  rename(ori = ORIGIN_STATE_NM) %>%
  rename(dest = DEST_STATE_NM) %>%
  ungroup()

#Need to remake since I actually have state names instead of abbreviations in this set
uniq_states <- union(t_DB1B_agg$ori,t_DB1B_agg$dest)
all_grid <- expand_grid(
  ori = uniq_states,
  dest = uniq_states,
  plus_one = 1
)

t_DB1B_agg<-right_join(t_DB1B_agg,all_grid,
                     by=join_by(ori,dest)) %>%
  mutate(pass = ifelse(is.na(pass),0,pass) + plus_one) %>% #Basically forces every pair to have 1 trip at minimum for RR calculation
  select(!plus_one)
t_DB1B_mirror <- t_DB1B_agg %>%
  rename(ori_mirror = dest) %>%
  rename(dest_mirror = ori) %>%
  rename(pass_mirror = pass)

df_rr_db1b <- left_join(t_DB1B_agg,t_DB1B_mirror,
                        by = join_by("ori" == "ori_mirror",
                                     "dest" == "dest_mirror")) %>%
  rename(x = ori) %>%
  rename(y = dest) %>%
  mutate(pass_xy = pass + pass_mirror,
         pass_total = sum(pass_xy)) %>%
  group_by(x) %>%
  mutate(pass_x = sum(pass_xy)) %>% ungroup() %>%
  group_by(y) %>%
  mutate(pass_y = sum(pass_xy)) %>% ungroup() %>%
  mutate(RR_air = pass_xy * pass_total / pass_x / pass_y) %>%
  select(x,y,RR_air)

write_csv(df_rr_db1b,"data/rr_air_db1b.csv")

combo_set <- df_rr_t100 %>%
  rename(RR_t100 = RR_air) %>%
  left_join(df_rr_db1b,
            by = join_by(x,y)) %>%
  rename(RR_db1b = RR_air) 

ggplot(combo_set) + 
  aes(x = RR_db1b,y=RR_t100) + 
  geom_point() + 
  geom_smooth() +
  xlim(0,10) +
  ylim(0,10)

#Quick summary of the datasets in comments: T100 has sine very extreme RR because it picks up small regional carriers doing commuter flights
#These are not a significant sample of the DB1B set so we find that DB1B data tends to be more in line with expected patterns
#Additioanally flights to and from hubs play a much larger role in the T100 data since it divides up a single trip into all of its legs
#For the purpose of identifying the travel patterns between states DB1B seems more helpful, although it will miss some of that within airport mixing
#But since sequencing will happen at the final destination, it probably wouldn't even catch that signal anyways