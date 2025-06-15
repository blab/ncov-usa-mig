library(tidyverse)

df_travel <- fread("data/interstate_travel_dot.csv")  %>%
  rename(x = Origin) %>%
  rename(y = Destination) %>% 
  rename(mode = "Pivot Field Names") %>%
  rename(trips_dir = "Pivot Field Values") %>% #Directional trips
  mutate(trips_dir = as.numeric(trips_dir)) #Switch from into to num since it will confuse the math later

df_travel$trips_xy <- 0 
STATE_LIST <- unique(df_travel$x)
for(i in STATE_LIST){
  for(j in STATE_LIST){
    if(i==j){
      df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir
    }else{
      df_travel[x == i & y == j]$trips_xy <- df_travel[x == i & y == j]$trips_dir + df_travel[x == j & y == i]$trips_dir
    }
  }
}
# Compute symmetric trips
df_travel<- df_travel %>%
  mutate(pair_id = pmap_chr(list(x, y), ~ paste(sort(c(..1, ..2)), collapse = "_"))) %>%
  group_by(pair_id) %>%
  mutate(trips_xy = sum(trips_dir, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-pair_id)

df_travel <- df_travel %>%
  group_by(x) %>%
  mutate(trips_x = sum(trips_xy)) %>%
  group_by(y) %>%
  mutate(trips_y = sum(trips_xy)) %>%
  ungroup() %>%
  mutate(trips_total = sum(trips_xy),
         same_state = (x == y), #Use to stratify for out of state travel only
         RR_trips = (trips_xy/trips_x)/(trips_y/trips_total))

write_csv(df_travel,"data/nhts_rr.csv")
# df_travel_out <- df_travel %>%
#   filter(!same_state) %>%
#   group_by(x) %>%
#   mutate(trips_x = sum(trips_xy)) %>%
#   group_by(y) %>%
#   mutate(trips_y = sum(trips_xy)) %>%
#   ungroup() %>%
#   mutate(trips_total = sum(trips_xy),
#          RR_trips = (trips_xy/trips_x)/(trips_y/trips_total)) %>%
#   rowwise() 
