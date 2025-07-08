#File: make_state_nb_list.R
#Author(s): Amin Bemanian
#Date: 8/19/24
#Description: R script that obtains state shapefiles and makes an neighbor (adjacency) and euclidean distance matrix using spdep/sf. 
#Does not need to be re-run with scenarios after initial run.
#TODO: Integrate with NHTS?

library(spdep)
library(tigris)
library(dplyr)
library(data.table)

options(tigris_use_cache = TRUE) #Save the shapefiles for future use
MAX_LAG <- 12 #Maximum higher order neighbors, California and Maine have the farthest queen adjacency distance at 11, may need to change if you use rook adjacency

#Okay to use generalized states file (cb=TRUE) since we just need this for an adjacency matrix
all_states <- states(cb=TRUE)
nb_states <- poly2nb(all_states,row.names = all_states$NAME) #Defaults to queen-adjacency
nb_lagged_states <- nblag(nb_states,MAX_LAG) 
STATE_NAMES <- all_states$NAME

#Switch from sparse lists to a long matrix

nb_matrix <- expand.grid(STATE_NAMES,STATE_NAMES,stringsAsFactors = FALSE) %>% as.data.table
colnames(nb_matrix) <- c("state_x","state_y")

nb_matrix$nb_dist <- as.numeric(-1) #Temporarily set default distance to -1 to initialize 

for(l in 1:MAX_LAG){
  nb_level <- nb_lagged_states[[l]]
  for(i in 1:length(STATE_NAMES)){
    name_i <- STATE_NAMES[i]
    nb_i <- nb_level[[i]] 
    #Set self to having a distance of 0, just do this in the first step
    if(l == 1){nb_matrix[state_x %in% name_i & state_y %in% name_i]$nb_dist <- 0}
    for(j in nb_i){
      name_j <- STATE_NAMES[j]
      nb_matrix[state_x %in% name_i & state_y %in% name_j]$nb_dist <- l
    }
  }
}

nb_matrix[nb_dist == -1]$nb_dist <- NA #Switch


#Second half is calculation of euclidean distances between state (centroids)
state_centroids <- st_centroid(all_states)
state_dists <- st_distance(state_centroids) * 1E-3 #Convert to km
attributes(state_dists)$units <- "km"
nb_matrix$euclid_dist <- -1.0 # Similarly set to -1 to initialize numeric values

for(i in 1:length(STATE_NAMES)){
  name_i <- STATE_NAMES[i]
  for(j in 1:length(STATE_NAMES)){
    name_j <- STATE_NAMES[j]
    nb_matrix[state_x %in% name_i & state_y %in% name_j]$euclid_dist <- state_dists[i,j]
  }
}


nb_matrix<-nb_matrix[order(nb_matrix$state_x,nb_matrix$state_y),]
readr::write_tsv(nb_matrix,"../data/nb_dist_states.tsv")

