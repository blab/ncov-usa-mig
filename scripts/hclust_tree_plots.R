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
library(RColorBrewer)
library(ggrepel)
library(ggdendro)
library(vegan)
library(usmap)


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
#Manual ordering for labels
DIV_ORDER<- c("Pacific","Mountain","West North Central","East North Central","Middle Atlantic",
              "New England","South Atlantic","East South Central","West South Central","Alaska","Hawaii","Pacific Territories","Caribbean Territories")

dir.create(paste0("figs/",scenario)) #Make a directory in case it doesn't already exist
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

STATE_CUT_HEIGHT <- 0.37

state_dist_mat <- exp(-state_rr_mat) %>% as.matrix() %>% as.dist()
state_hc <- hclust(state_dist_mat)
state_hclust <- cutree(state_hc,h=STATE_CUT_HEIGHT) %>% 
  as_tibble(rownames = "state") %>%
  rename(cluster = value) %>%
  mutate(cluster=factor(cluster))

state_dendro <- dendro_data(state_hc)
state_dendro_stems <- state_dendro$segments
state_dendro_leaves <- state_dendro$labels %>%
  rename(state = label) %>%
  inner_join(state_hclust,by="state")

fig_state_dendro <- ggplot() +
  theme_void() +
  geom_segment(data=state_dendro_stems,
               aes(x=x,y=y,xend=xend,yend=yend)) +
  geom_point(data=state_dendro_leaves,
             aes(x=x,y=y,color=cluster)) +
  geom_text(data=state_dendro_leaves,
            aes(x=x,y=y,label=state,color=cluster),
            hjust = -0.1,
            show.legend = FALSE) +
  geom_hline(yintercept = STATE_CUT_HEIGHT,linetype=2) +
  labs(title = "State Dendrogram",
       color= "Cluster Number") +
  scale_y_reverse() +
  coord_flip() +
  ylim(0.6,-0.2) +
  theme(plot.title = element_text(hjust = 0.5))

fn_state_dedro <- paste0("figs/",scenario,"/state_dendrogram.jpg")
ggsave(fn_state_dedro,plot=fig_state_dendro,width=2400,height=2400,units="px")

set.seed(17)
state_NMDS <- metaMDS(state_dist_mat, k=2)
state_PCoA <- cmdscale(state_dist_mat, eig = TRUE)
state_PCoA$perc_explained <- state_PCoA$eig/sum(state_PCoA$eig) * 100 #Add a percent explained

state_kclust_nmds <- kmeans(scores(state_NMDS),centers = 10)$cluster %>%
  as_tibble(rownames="state") %>%
  rename(cluster = value) %>%
  mutate(cluster = as.factor(cluster))

state_kclust_pcoa <- kmeans(state_PCoA$points,centers = 10)$cluster %>%
  as_tibble(rownames="state") %>%
  rename(cluster = value) %>%
  mutate(cluster = as.factor(cluster))

theme_set(theme_bw())

fig_state_NMDS_hclust <- scores(state_NMDS) %>%
  as_tibble(rownames="state") %>%
  inner_join(state_hclust, by="state") %>% 
  ggplot(aes(x=NMDS1,y=NMDS2)) +
  geom_point(aes(color=cluster)) +
  geom_text_repel(aes(label=state,color=cluster),max.overlaps = 15,
                  show.legend = FALSE) +
  labs(title = "Non-metric Multidimensional Scaling of States",
       x = "Component #1",
       y = "Component #2",
       color = "Hierarchical Cluster") +
  theme(plot.title = element_text(hjust = 0.5))
fn_state_NMDS_hclust <- paste0("figs/",scenario,"/state_NMDS_hclust.jpg")
ggsave(fn_state_NMDS_hclust,plot=fig_state_NMDS_hclust,width=3000,height=2400,units="px")

fig_state_NMDS_kclust <- scores(state_NMDS) %>%
  as_tibble(rownames="state") %>%
  inner_join(state_kclust_nmds, by="state") %>% 
  ggplot(aes(x=NMDS1,y=NMDS2)) +
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

fig_state_PCoA_hclust <- scores(state_PCoA) %>%
  as_tibble(rownames="state") %>%
  inner_join(state_hclust, by="state") %>% 
  ggplot(aes(x=Dim1,y=Dim2)) +
  geom_point(aes(color=cluster))  +
  geom_text_repel(aes(label=state,color=cluster),max.overlaps = 20,
                  show.legend = FALSE) +
  labs(title = "Principal Coordinate Analysis of States",
       x = paste0("Component #1 Percent Explained: ",round(state_PCoA$perc_explained[1],1),"%"),
       y = paste0("Component #2 Percent Explained: ",round(state_PCoA$perc_explained[2],1),"%"),
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5))
fn_state_PCoA_hclust <- paste0("figs/",scenario,"/state_PCoA_hclust.jpg")
ggsave(fn_state_PCoA_hclust,plot=fig_state_PCoA_hclust,width=3000,height=2400,units="px")

fig_state_PCoA_kclust <- scores(state_PCoA) %>%
  as_tibble(rownames="state") %>%
  inner_join(state_kclust_pcoa,by="state") %>% 
  ggplot(aes(x=Dim1,y=Dim2)) +
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


map_state_hclust <- plot_usmap(regions="states",data=state_hclust,values="cluster") +
  labs(title = "State Map of Hierarchical Clusters",
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5))
fn_map_state_hclust <- paste0("figs/",scenario,"/maps/state_hclust.jpg")
ggsave(fn_map_state_hclust,plot=map_state_hclust,width=2400,height=1800,units="px")

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

df_clust_state <- inner_join(state_hclust,state_kclust_nmds,by=join_by(state)) %>%
  rename(hclust = cluster.x) %>%
  rename(kclust_nmds = cluster.y) %>%
  inner_join(state_kclust_pcoa, by=join_by(state)) %>%
  rename(kclust_pcoa = cluster)

fn_clust_state <- paste0("results/",scenario,"/df_state_clusters.tsv")
readr::write_tsv(clust_state,fn_clust_state)
