#SIR Patch Model Simulation - Single Strain
#Will simulate transmission of a new virus using a potch model
#Will assume patches are on a 21x21 grid, with 1000 individuals each
#Patch number will assign to coordiantes: i = x + 21*y
#i.e. x = i mod 21 and y = floor(21)
#Variable of interest is proportion of population in I
#set.seed(17)
library(ggplot2)
library(tidyverse)
library(spdep)
library(magick)

#Naming for plot file names and  titles
OUT_FN <- "multi"
OUT_TITLE <- "Multiple Sources"
SAVE_FIGS <- FALSE

#Spatial parameters
NUM_ROW <- 21
NUM_COL <- 21
NB_WT <- 0.25
JUMP_WT <- NB_WT/5
num_patches <- NUM_ROW * NUM_COL

coord_to_index <- function(x,y){
  i <- x + NUM_ROW * y + 1
  return(i)
}

index_to_coord <- function(i){
  x <- (i-1) %% NUM_COL
  y <- floor((i-1)/NUM_ROW)
  return(c(x,y))
}

SPECIAL_WT <- list(c(11,252),c(240,426), c(106,368),c(67,423))

calc_nb_wts <- function(i,incl_special=TRUE){
  coord <- index_to_coord(i)
  x <- coord[1]
  y <- coord[2]
  wt <- rep(0,num_patches)
  wt[i] <- 1
  if(x > 0){wt[coord_to_index(x-1,y)]<-NB_WT}
  if(x < (NUM_COL-1)){wt[coord_to_index(x+1,y)]<-NB_WT}
  if(y > 0){wt[coord_to_index(x,y-1)]<-NB_WT}
  if(y < (NUM_ROW-1)){wt[coord_to_index(x,y+1)]<-NB_WT}
  if(length(SPECIAL_WT) > 0 & incl_special){for(p in 1:length(SPECIAL_WT)){
    pair <-SPECIAL_WT[[p]]
    if(i == pair[1]){wt[pair[2]]<-NB_WT}
    if(i == pair[2]){wt[pair[1]]<-NB_WT}
  }}
  wt<-wt/sum(wt)
  wt
}

wt_matrix <- do.call(rbind,lapply(1:num_patches,calc_nb_wts))
w <- wt_matrix
for(i in 1:num_patches){
  w[i,i] <- 0
  w[i,] <- w[i,]/sum(w[i,])
}
lw<-mat2listw(w, style = "W")

simple_wt_matrix <- do.call(rbind,lapply(1:num_patches,calc_nb_wts,incl_special=FALSE))
w_s <- simple_wt_matrix
for(i in 1:num_patches){
  w_s[i,i] <- 0
  w_s[i,] <- w_s[i,]/sum(w_s[i,])
}
lw_s<-mat2listw(w_s, style = "W")

#Sim parameters
#Values of these are picked to provide enough time
#For infection wave to propagate and be visualized 
R0 <- 10  
RECOVERY_TIME <- 10
T_MAX <- 150
BETA_SD <- 0.5
SEED_INDICES <- runif(1,min=1,max=num_patches) %>% round() #Starting indices
#Assume patch population size and starting infection size is fixed
TOTAL_POP <- 1000 
STARTING_INF <- 20

S_start <- rep(TOTAL_POP,num_patches)
S_start[SEED_INDICES] <- TOTAL_POP - STARTING_INF
I_start <- rep(0,num_patches)
I_start[SEED_INDICES] <- STARTING_INF
R_start <- rep(0,num_patches) #Totally naive population

ARRIVAL_TRHESHOLD <- 0.1
arrival_vector <- rep(T_MAX+1,num_patches) #Vector of arrival time points

#Dependent constants
gamma_inf <- 1/RECOVERY_TIME
beta_inf <- R0 * gamma_inf
N <- S_start + I_start + R_start
beta_scale <- exp(rnorm(n=num_patches,mean = 0, sd = BETA_SD)) #Adds patch variance to beta

#For the rate equations of each patch we assume no one permanently relocates
#Therefore the susceptible pool will just be the susceptible individuals in the patch
#The infectious pool will be the weighted sum of the infectious pools that come into contact with the patch
#For recovery however, it will again only be the infectious individuals residing in the patch 

delta_s <- function(s_patch,i_patch,n_patch,beta_inf){
  as.numeric(-beta_inf*beta_scale/n_patch * i_patch %*% wt_matrix * s_patch)
}
delta_i <- function(s_patch,i_patch,n_patch,beta_inf,gamma_inf){
  as.numeric(beta_inf*beta_scale/n_patch * i_patch %*% wt_matrix * s_patch  - gamma_inf * i_patch)
}
delta_r <- function(i_patch,gamma_inf){
  i_patch * gamma_inf
}

d <- tibble(t = 0, 
           s = list(S_start),
           i = list(I_start),
           r = list(R_start)
           )

for(t in 1:T_MAX){
  s <- (d$s[[t]] + delta_s(d$s[[t]],d$i[[t]],N,beta_inf)) %>% sapply(max,0) %>% sapply(min,N)
  i <- (d$i[[t]] + delta_i(d$s[[t]],d$i[[t]],N,beta_inf,gamma_inf)) %>% sapply(max,0) %>% sapply(min,N)
  r <- (d$r[[t]] + delta_r(d$i[[t]],gamma_inf)) %>% sapply(max,0) %>% sapply(min,N)
  x <- tibble(t,s=list(s),i=list(i),r=list(r))
  d <- add_row(d,x)
  
  above_threshold <- which(i > ARRIVAL_TRHESHOLD*N)
  arrival_vector[above_threshold] <- arrival_vector[above_threshold] %>% sapply(min,t)
}
arrival_vector<- arrival_vector %>% unlist() #Unsure why this is making a list

plot_prop <- function(p,title){
  coords <- t(sapply(1:num_patches,index_to_coord))
  df <- tibble(
    x = coords[,1],
    y = coords[,2],
    val = p
  )
  ggplot(df, aes(x = x, y = y, fill = p)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "yellow", high = "red",limits=c(0,1)) +
#   geom_text(aes(label=round(p,3)),size=3) +
    coord_fixed() +
    scale_y_reverse() +  # optional: flip vertically to make (0,0) bottom-left
    theme_minimal() +
    labs(title = title, fill = "Proportion")
}

#Adapted from nagraj.net
make_gif<-function(frame_out,file_out,FPS=25){
  imgs <- list.files(frame_out, full.names = TRUE)
  img_list <- imgs %>% lapply(function(i){
    image_read(i)
  })
  img_joined <- image_join(img_list)
  img_animated <- image_animate(img_joined, fps = FPS,optimize=TRUE)
  image_write(image=img_animated,path=file_out,quality=100)
}

plot_arrival <- function(t,title){
  coords <- t(sapply(1:num_patches,index_to_coord))
  df <- tibble(
    x = coords[,1],
    y = coords[,2],
    val = t
  )
  ggplot(df, aes(x = x, y = y, fill = t)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "khaki", high = "dodgerblue") +
    #geom_text(aes(label=round(t,3)),size=3) +
    coord_fixed() +
    scale_y_reverse() +  # optional: flip vertically to make (0,0) bottom-left
    theme_minimal() +
    labs(title = title, fill = "Time")
}


moran_bounds <- moran.mc(d$i[[T_MAX/2]]+rnorm(num_patches,sd=100),
                         lw,
                         alternative = "two.sided",
                         nsim=1E5)
moran_ub <- quantile(moran_bounds$res,0.975)
moran_lb <- quantile(moran_bounds$res,0.025)

inf_moran <- function(t){
  i<-d$i[[t]]
  err <- rnorm(num_patches,sd=10) #Simulate measurement error
  m <- moran.mc(i+err,lw,nsim=1)$statistic
  tibble(t,m)
}
df_moran<-map_dfr(1:T_MAX,inf_moran) 

moran_arrival <- moran.mc(round(arrival_vector+rnorm(num_patches,0,2)),
                          lw,
                          alternative = "two.sided",
                          nsim = 1E4)

if(SAVE_FIGS){
  for(t in 1:T_MAX){
    plot_prop(d$r[[t]]/N,paste0("Recovered at time: ",t))
    ggsave(sprintf("figs/sim/frames/recovered/frame_%03d.jpg", t),
           height = 6,
           width = 6,
           units = "in",
           dpi = 150)
  }
  
  
  for(t in 1:T_MAX){
    plot_prop(d$i[[t]]/N,paste0("Infected at time: ",t))
    ggsave(sprintf("figs/sim/frames/infected/frame_%03d.jpg", t),
           height = 6,
           width = 6,
           units = "in",
           dpi = 150)
  }
  
  make_gif("figs/sim/frames/recovered",
           paste0("figs/sim/recovered_simulation_",OUT_FN,".gif"))
  make_gif("figs/sim/frames/infected",
           paste0("figs/sim/infected_simulation_",OUT_FN,".gif"))
  
  ggplot(df_moran,aes(x=t,y=m)) + 
    geom_line() + ggplot2::ylim(-0.1,1) +
    geom_ribbon(aes(ymin=moran_lb,ymax=moran_ub),alpha=0.2) +
    theme_bw() + 
    labs(title = OUT_TITLE)
  ggsave(paste0("figs/sim/moran_",OUT_FN,".jpg"),
         width = 8,
         height = 5,
         units = "in",
         dpi = 150)
  
  plot_arrival(arrival_vector,paste0("Arrival Time Heatmap - ",OUT_TITLE)) +
    labs(subtitle = paste0("Moran I: ", round(moran_arrival$statistic,3), 
                           " - p-val: ", round(moran_arrival$p.value,3)))
  ggsave(paste0("figs/sim/arrival_",OUT_FN,".jpg"),
         width = 6,
         height = 6,
         units = "in",
         dpi = 150
  )
}

arrival_title <- paste0("Origin Infection: (",index_to_coord(SEED_INDICES)[1],", ",index_to_coord(SEED_INDICES)[2],")")
plot_arrival(arrival_vector,paste0("Arrival Time Heatmap - ",arrival_title)) +
  labs(subtitle = paste0("Moran I: ", round(moran_arrival$statistic,3), 
                         " - p-val: ", round(moran_arrival$p.value,3)))

stop("End of main script")

#Code to look at interoperator reliability
clust69 <- hclust(dist(df_moving_arrivals$ori69)) %>% cutree(k=4)
clust92 <- hclust(dist(df_moving_arrivals$ori92)) %>% cutree(k=4)
clust198 <- hclust(dist(df_moving_arrivals$ori198)) %>% cutree(k=4)
clust309 <- hclust(dist(df_moving_arrivals$ori309)) %>% cutree(k=4)
mclust::adjustedRandIndex(clust69,clust92)
mclust::adjustedRandIndex(clust69,clust198)
mclust::adjustedRandIndex(clust92,clust198)
mclust::adjustedRandIndex(clust69,clust309)
mclust::adjustedRandIndex(clust92,clust309)
mclust::adjustedRandIndex(clust198,clust309)


#Code for making example lagged plots for salon
test_prop <- ((round(d$i[[30]]) * rnorm(num_patches,1,0.3)) %>% sapply(min,1E3))/1E3
plot_prop(test_prop,title="Infected")
lagged_prop <- w_s %*% test_prop
plot_prop(lagged_prop,title="Lagged Infected")
df_test <- tibble(prop=test_prop,lagged_prop=w_s %*% test_prop)
ggplot(df_test,aes(x=test_prop,y=lagged_prop)) + 
  geom_point() + 
  geom_smooth(method=lm) +
  labs(x="Infected Proportion",y="Lagged Proportion") +
  theme_bw()

lm(lagged_prop ~ test_prop,data=df_test)
