library(sf)
library(tidyverse)
library(patchwork)

df_RR_state <- read_tsv("results/CAM_1000/df_RR_by_state.tsv")
cam_map <- st_read("data/shp-files/CAM-simplified.gpkg")
st_bbox(cam_map)
st_crs(cam_map)
cam_map <- cam_map %>% 
  st_transform(st_crs("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-100 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
st_crs(cam_map)
st_bbox(cam_map)

ORI <- "Illinois"

#Set as bounds for RR and transform to log for display purposes
UB<-2 
LB <- 0.5 
fill_bound <- function(x){
  log_x <- log10(x)
  max(min(log_x,log10(UB)),log10(LB)) 
}

#Subset to state to plot
df_subset <- df_RR_state %>% 
  filter(x == ORI) %>%
  rename(state = y)

df_map <- cam_map %>% 
  left_join(df_subset,by = join_by(NAME_En == state))
df_map$logRR <- sapply(df_map$RR,fill_bound)

df_main <- df_map %>% filter(NAME_En != "Hawaii") #Mainland
df_hi <- df_map %>% filter(NAME_En == "Hawaii") #Hawai'i

fill_scale <- scale_fill_gradient2(
  name = "RR",
  high = "#D67C34",
  low = "#4C90C0",
  limits = c(log10(LB), log10(UB)),
  breaks = c(log10(LB), log10((1+LB)/2), log10(1), log10((1+UB)/2), log10(UB)),
  labels = c(LB, (1+LB)/2, 1, (1+UB)/2, UB)
)

p_main <- ggplot(df_main) + 
  geom_sf(aes(fill = logRR), color = "black", size = 1) +
  fill_scale +
  theme_minimal() +
  labs(title = paste0("Enrichment RRs for ", ORI))

p_hi <- ggplot(df_hi) + 
  geom_sf(aes(fill = logRR), color = "black", size = 1) +
  fill_scale +
  theme_void()

ggplot(df_map) + 
  geom_sf(aes(fill=logRR),color="black",size =1) +
  scale_fill_gradient2(name="RR",
                       high = "#D67C34",
                       low = "#4C90C0",
                       limits=c(log10(LB),log10(UB)),
                       breaks =c(log10(LB),log10((1+LB)/2),log10(1),log10((1+UB)/2),log10(UB)),
                       labels = c(LB,(1+LB)/2,1,(1+UB)/2,UB)) +
  theme_minimal() +
  labs(title = paste0("Enrichment RRs for ",ORI))
ggsave("figs/test_map.jpg",height = 2000,width=1600,units = "px",dpi=192)
