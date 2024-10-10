#REQUIRES RESTANDARDIZATION WITH THE OTHER SCRIPT FILES

library(usmap)
library(dplyr)
library(ggplot2)

load("./data/df_variant_div.Rdata")

var_to_map <- "XBB.1.5"
week_to_map <- lubridate::floor_date(as.Date("2022-12-01"),"week") #Rounds to the Sunday of that week so it matches with the week bins
state_perc <-  df_variant_div %>%
    filter(week == week_to_map) %>%
    filter(variant == var_to_map) %>%
    mutate(state = division) %>% #To work with USMaps package
    select(state,perc_voi)
head(state_perc)

var_usmap <- plot_usmap(data=state_perc,
    values = "perc_voi") +
    labs(title = paste("Variant:",var_to_map,"Week:",week_to_map)) +
    scale_fill_continuous(limits=c(0,100),
                          name="Variant Frequency",
                          type="gradient",
                          low="lightgoldenrod",
                          high="red3") +
    theme(legend.position=c(0.9,0.05))

bitmap(paste("./maps/Var_USMap_",var_to_map,"_",week_to_map,".png",sep=""),
       res=600,height=4,width=7)
print(var_usmap)
dev.off()
