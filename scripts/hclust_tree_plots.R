#File: state_analysis.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Make geographical hierarchical cluster plots
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)
library(ggdendro)
library(vegan)
library(usmap)
library(cluster)
library(tsne)
library(factoextra) #Needed for elbow plots


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--dist', type = 'character', help = "Cluster Distance Function Reciprocal (R) or Neg-Expontential (E)", default = "E")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
dist_type <- args$dist #'E' or 'R'

#Manual ordering for labels
DIV_ORDER<- c("Pacific","Mountain","West North Central","East North Central","Middle Atlantic",
              "New England","South Atlantic","East South Central","West South Central","Alaska","Hawaii","Pacific Territories","Caribbean Territories")

dir.create(paste0("figs/",scenario,"/clust")) #Make a directory in case it doesn't already exist
fn_div_rr <- paste0("results/",scenario,"/df_RR_by_census_div.tsv")
div_rr <- fread(fn_div_rr)

div_rr_mat <- div_rr %>% select(c("x","y","RR")) %>%
  reshape(idvar="x",v.names = c("RR"),timevar = "y",direction = "wide") 
div_names <- div_rr_mat$x
div_rr_mat <- div_rr_mat[,2:ncol(div_rr_mat)]
rownames(div_rr_mat) <- div_names
colnames(div_rr_mat) <- div_names

div_dist_mat <- exp(-div_rr_mat) %>% as.matrix()
rownames(div_dist_mat) <- div_names
div_dist_mat[is.na(div_dist_mat)] <- 1E1
div_dist_mat <- as.dist(div_dist_mat)
div_hc <- hclust(div_dist_mat)

fn_state_rr <- paste0("results/",scenario,"/df_RR_by_state.tsv")
state_rr <- fread(fn_state_rr)

state_rr_mat <- state_rr %>% select(c("x","y","RR")) %>%
  reshape(idvar="x",v.names = c("RR"),timevar = "y",direction = "wide") 
state_names <- state_rr_mat$x
state_rr_mat <- state_rr_mat[,2:ncol(state_rr_mat)]
rownames(state_rr_mat) <- state_names
colnames(state_rr_mat) <- state_names

#STATE_CUT_HEIGHT <- 0.37
NUM_CLUSTS <- 6
METHOD_LIST <- c("single","complete","ward.D2","average","centroid")

state_dist_mat <- case_when(
  dist_type == 'E' ~ exp(-state_rr_mat) %>% as.matrix() %>% as.dist(),
  dist_type == 'R' ~ 1/(state_rr_mat) %>% as.matrix() %>% as.dist(),
  TRUE ~ NA
)

dist_full_name <- case_when( #For title purposes
  dist_type == 'E' ~ "Negative Exponential",
  dist_type == 'R' ~ "Reciprocal",
  TRUE ~ NA
)


for(m in METHOD_LIST){
  
  subtitle_text <- paste("Agglomerative",m,"-",dist_full_name,"Distance") 
  state_hc <- hclust(state_dist_mat,method=m)
  state_hclust <- cutree(state_hc,k=NUM_CLUSTS) %>% 
    as_tibble(rownames = "state") %>%
    rename(cluster = value) %>%
    mutate(cluster=factor(cluster))

  state_dendro <- dendro_data(state_hc)
  state_dendro_stems <- state_dendro$segments
  state_dendro_leaves <- state_dendro$labels %>%
    rename(state = label) %>%
    inner_join(state_hclust,by="state")
  

  fig_state_dendro <- ggplot() +
    theme_bw() +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    geom_segment(data=state_dendro_stems,
                aes(x=x,y=y,xend=xend,yend=yend)) +
    geom_point(data=state_dendro_leaves,
              aes(x=x,y=y,color=cluster)) +
    geom_text(data=state_dendro_leaves,
              aes(x=x,y=y,label=state,color=cluster),
              hjust = -0.1,
              show.legend = FALSE) +
#    geom_hline(yintercept = STATE_CUT_HEIGHT,linetype=2) 
    labs(title = "State Dendrogram",
        subtitle = subtitle_text,
        color = "Cluster Number",
        y = "Cluster Distance") +
    scale_y_reverse() +
    coord_flip() +
    ylim(max(state_dendro_stems$y),-0.25*max(state_dendro_stems$y)) + #Dynamic range that lets labels fit
    theme(plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "bottom")

  #fn_state_dedro <- paste0("figs/",scenario,"/state_dendrogram_"dist_type,"_agg_",m,".jpg")
  #ggsave(fn_state_dedro,plot=fig_state_dendro,width=2400,height=2400,units="px")

  map_state_hclust <- plot_usmap(regions="states",data=state_hclust,values="cluster") +
  labs(title ="State Map of Hierarchical Clusters",
       subtitle = subtitle_text,
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom")

  #fn_map_state_hclust <- paste0("figs/",scenario,"/maps/state_hclust_"dist_type,"_agg_",m,".jpg")
  #ggsave(fn_map_state_hclust,plot=map_state_hclust,width=2400,height=1800,units="px")
  
  #Combined dendogram and map
  fn_state_combo_hclust <- paste0("figs/",scenario,"/clust/state_hclust_",dist_type,"_agg_",m,".jpg")
  fig_state_combo_hclust <- ggarrange(fig_state_dendro,map_state_hclust,nrow=2)
  ggsave(fn_state_combo_hclust, fig_state_combo_hclust,width=2400,height=4800,units="px")
}

#Iterate over number of clusters
for(k_clust in 2:8){
  m <- "ward.D2"
  subtitle_text <- paste("Agglomerative",m,"-",dist_full_name,"Distance") 
  state_hc <- hclust(state_dist_mat,method=m)
  state_hclust <- cutree(state_hc,k=k_clust) %>% 
    as_tibble(rownames = "state") %>%
    rename(cluster = value) %>%
    mutate(cluster=factor(cluster))
  
  state_dendro <- dendro_data(state_hc)
  state_dendro_stems <- state_dendro$segments
  state_dendro_leaves <- state_dendro$labels %>%
    rename(state = label) %>%
    inner_join(state_hclust,by="state")
  
  
  fig_state_dendro <- ggplot() +
    theme_bw() +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    geom_segment(data=state_dendro_stems,
                 aes(x=x,y=y,xend=xend,yend=yend)) +
    geom_point(data=state_dendro_leaves,
               aes(x=x,y=y,color=cluster)) +
    geom_text(data=state_dendro_leaves,
              aes(x=x,y=y,label=state,color=cluster),
              hjust = -0.1,
              show.legend = FALSE) +
    #    geom_hline(yintercept = STATE_CUT_HEIGHT,linetype=2) 
    labs(title = "State Dendrogram",
         subtitle = subtitle_text,
         color = "Cluster Number",
         y = "Cluster Distance") +
    scale_y_reverse() +
    coord_flip() +
    ylim(max(state_dendro_stems$y),-0.25*max(state_dendro_stems$y)) + #Dynamic range that lets labels fit
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  #fn_state_dedro <- paste0("figs/",scenario,"/state_dendrogram_"dist_type,"_agg_",m,".jpg")
  #ggsave(fn_state_dedro,plot=fig_state_dendro,width=2400,height=2400,units="px")
  
  map_state_hclust <- plot_usmap(regions="states",data=state_hclust,values="cluster") +
    labs(title ="State Map of Hierarchical Clusters",
         subtitle = subtitle_text,
         color = "Cluster Number") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  #fn_map_state_hclust <- paste0("figs/",scenario,"/maps/state_hclust_"dist_type,"_agg_",m,".jpg")
  #ggsave(fn_map_state_hclust,plot=map_state_hclust,width=2400,height=1800,units="px")
  
  #Combined dendogram and map
  fn_state_combo_hclust <- paste0("figs/",scenario,"/clust/state_hclust_",k_clust,"_clusts_",m,".jpg")
  fig_state_combo_hclust <- ggarrange(fig_state_dendro,map_state_hclust,nrow=2)
  ggsave(fn_state_combo_hclust, fig_state_combo_hclust,width=2400,height=4800,units="px")
}

fviz_nbclust(
  exp(-state_rr_mat),
  hcut,
  method='wss',
  k.max = 8
)

set.seed(17)
state_NMDS <- metaMDS(state_dist_mat, k=2)
state_PCoA <- cmdscale(state_dist_mat, eig = TRUE)
state_PCoA$perc_explained <- state_PCoA$eig/sum(state_PCoA$eig) * 100 #Add a percent explained
state_tSNE <- tsne(state_dist_mat,k=2) %>% data.frame() %>% magrittr::set_rownames(state_names)

fviz_nbclust(scores(state_NMDS),
             kmeans,
             method='silhouette',
             k.max = 8
)

state_kclust_nmds <- kmeans(scores(state_NMDS),centers = 3)$cluster %>%
  as_tibble(rownames="state") %>%
  rename(cluster = value) %>%
  mutate(cluster = as.factor(cluster))

fviz_nbclust(state_PCoA$points,
             kmeans,
             method='silhouette',
             k.max = 8
)


state_kclust_pcoa <- kmeans(state_PCoA$points,centers = 3)$cluster %>%
  as_tibble(rownames="state") %>%
  rename(cluster = value) %>%
  mutate(cluster = as.factor(cluster))

fviz_nbclust(state_tSNE,
             kmeans,
             method='silhouette',
             k.max = 8
)

state_kclust_tSNE <- kmeans(state_tSNE,centers=3)$cluster %>%
  as_tibble(rownames="state") %>%
  rename(cluster = value) %>%
  mutate(cluster = as.factor(cluster))

print("Completed clustering")

theme_set(theme_bw())

# fig_state_NMDS_hclust <- scores(state_NMDS) %>%
#   as_tibble(rownames="state") %>%
#   inner_join(state_hclust, by="state") %>% 
#   ggplot(aes(x=NMDS1,y=NMDS2)) +
#   geom_point(aes(color=cluster)) +
#   geom_text_repel(aes(label=state,color=cluster),max.overlaps = 15,
#                   show.legend = FALSE) +
#   labs(title = "Non-metric Multidimensional Scaling of States",
#        x = "Component #1",
#        y = "Component #2",
#        color = "Hierarchical Cluster") +
#   theme(plot.title = element_text(hjust = 0.5))
# fn_state_NMDS_hclust <- paste0("figs/",scenario,"/state_NMDS_hclust.jpg")
# ggsave(fn_state_NMDS_hclust,plot=fig_state_NMDS_hclust,width=3000,height=2400,units="px")

fig_state_NMDS_kclust <- scores(state_NMDS) %>%
  as_tibble(rownames="state") %>%
  inner_join(state_kclust_nmds, by="state") %>% 
  ggplot(aes(x=-NMDS1,y=-NMDS2)) +
  geom_point(aes(color=cluster)) +
  geom_text_repel(aes(label=state,color=cluster),max.overlaps = 15,
                  show.legend = FALSE) +
  labs(title = "Non-metric Multidimensional Scaling of States",
       x = "Component #1",
       y = "Component #2",
       color = "K-means Cluster") +
  theme(plot.title = element_text(hjust = 0.5))
fn_state_NMDS_kclust <- paste0("figs/",scenario,"/state_NMDS_kclust.jpg")
ggsave(fn_state_NMDS_kclust,plot=fig_state_NMDS_kclust,width=3000,height=2400,units="px")

# fig_state_PCoA_hclust <- scores(state_PCoA) %>%
#   as_tibble(rownames="state") %>%
#   inner_join(state_hclust, by="state") %>% 
#   ggplot(aes(x=Dim1,y=Dim2)) +
#   geom_point(aes(color=cluster))  +
#   geom_text_repel(aes(label=state,color=cluster),max.overlaps = 20,
#                   show.legend = FALSE) +
#   labs(title = "Principal Coordinate Analysis of States",
#        x = paste0("Component #1 Percent Explained: ",round(state_PCoA$perc_explained[1],1),"%"),
#        y = paste0("Component #2 Percent Explained: ",round(state_PCoA$perc_explained[2],1),"%"),
#        color = "Cluster Number") +
#   theme(plot.title = element_text(hjust = 0.5))
# fn_state_PCoA_hclust <- paste0("figs/",scenario,"/state_PCoA_hclust.jpg")
# ggsave(fn_state_PCoA_hclust,plot=fig_state_PCoA_hclust,width=3000,height=2400,units="px")

fig_state_PCoA_kclust <- scores(state_PCoA) %>%
  as_tibble(rownames="state") %>%
  inner_join(state_kclust_pcoa,by="state") %>% 
  ggplot(aes(x=Dim1,y=-Dim2)) +
  geom_point(aes(color=cluster))  +
  geom_text_repel(aes(label=state,color=cluster),max.overlaps = 20,
                  show.legend = FALSE) +
  labs(title = "Principal Coordinate Analysis of States",
       x = paste0("Component #1 Percent Explained: ",round(state_PCoA$perc_explained[1],1),"%"),
       y = paste0("Component #2 Percent Explained: ",round(state_PCoA$perc_explained[2],1),"%"),
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5))
fn_state_PCoA_kclust <- paste0("figs/",scenario,"/state_PCoA_kclust.jpg")
ggsave(fn_state_PCoA_kclust,plot=fig_state_PCoA_kclust,width=3000,height=2400,units="px")

fig_state_tSNE_kclust <- scores(state_tSNE) %>%
  as_tibble(rownames="state") %>%
  inner_join(state_kclust_tSNE,by="state") %>% 
  ggplot(aes(x=-X1,y=-X2)) +
  geom_point(aes(color=cluster))  +
  geom_text_repel(aes(label=state,color=cluster),max.overlaps = 20,
                  show.legend = FALSE) +
  labs(title = "t-SNE Clustering of States",
       x = paste0("Dimension #1"),
       y = paste0("Dimension #2"),
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5))
fn_state_tSNE_kclust <- paste0("figs/",scenario,"/state_tSNE_kclust.jpg")
ggsave(fn_state_tSNE_kclust,plot=fig_state_tSNE_kclust,width=3000,height=2400,units="px")

map_state_PCoA_kclust <- plot_usmap(regions="states",data=state_kclust_pcoa,values="cluster") +
  labs(title = "State Map of PCoA Clusters",
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5))
fn_map_state_PCoA_kclust <- paste0("figs/",scenario,"/maps/state_PCoA_kclust.jpg")
ggsave(fn_map_state_PCoA_kclust,plot=map_state_PCoA_kclust,width=2400,height=1800,units="px")

map_state_NMDS_kclust <- plot_usmap(regions="states",data=state_kclust_nmds,values="cluster") +
  labs(title = "State Map of NMDS Clusters",
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5))
fn_map_state_NMDS_kclust <- paste0("figs/",scenario,"/maps/state_NMDS_kclust.jpg")
ggsave(fn_map_state_NMDS_kclust,plot=map_state_NMDS_kclust,width=2400,height=1800,units="px")

map_state_tSNE_kclust <- plot_usmap(regions="states",data=state_kclust_tSNE,values="cluster") +
  labs(title = "State Map of t-SNE Clusters",
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5))
fn_map_state_tSNE_kclust <- paste0("figs/",scenario,"/maps/state_tSNE_kclust.jpg")
ggsave(fn_map_state_NMDS_kclust,plot=map_state_NMDS_kclust,width=2400,height=1800,units="px")

# df_clust_state <- inner_join(state_hclust,state_kclust_nmds,by=join_by(state)) %>%
#   rename(hclust = cluster.x) %>%
#   rename(kclust_nmds = cluster.y) %>%
#   inner_join(state_kclust_pcoa, by=join_by(state)) %>%
#   rename(kclust_pcoa = cluster)

# fn_clust_state <- paste0("results/",scenario,"/df_state_clusters.tsv")
# readr::write_tsv(clust_state,fn_clust_state)
