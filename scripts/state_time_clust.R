#File: state_time_clust.R
#Author(s): Amin Bemanian
#Date: 8/13/24
#Description: Make geographical hierarchical cluster plots
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(tidyverse)
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
library(factoextra) 
library(jsonlite)
library(gridGraphics)


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--dist', type = 'character', help = "Cluster Distance Function Reciprocal (R) or Neg-Expontential (E)", default = "E")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- "USA_mini" #args$scenario
dist_type <- args$dist #'E' or 'R'


fn_path_rr <- paste0("results/",scenario,"/time_state/")
state_rr_all <- fread(paste0(fn_path_rr,"df_state_rr_all.tsv"))

state_rr_all_mat <- state_rr_all %>% select(c("x","y","RR")) %>%
  rename(name = x) %>%
  pivot_wider(names_from = y, values_from = RR) %>%
  arrange(name) %>%
  select(name, sort(setdiff(names(.), "name"))) %>%
  column_to_rownames("name") %>%
  as.matrix()

state_dist_mat <- case_when(
  dist_type == 'E' ~ exp(-state_rr_all_mat) %>% as.matrix() %>% as.dist(),
  dist_type == 'R' ~ 1/(state_rr_all_mat) %>% as.matrix() %>% as.dist(),
  TRUE ~ NA
)

dist_full_name <- case_when( #For title purposes
  dist_type == 'E' ~ "Negative Exponential",
  dist_type == 'R' ~ "Reciprocal",
  TRUE ~ NA
)

#Figuring out cluster height to cut at based on the full AGNES tree
hc_agnes_full <- agnes(state_dist_mat,diss = TRUE) %>% as.hclust()
h_a <- hc_agnes_full$height
t <- 1.25
cut_a <- mean(h_a) + t * sd(h_a)
cutree(hc_agnes_full, h = cut_a)

#Use silhouette plots to confirm reasonable height on the full set
fviz_nbclust(x=state_rr_all_mat,diss=state_dist_mat,FUN=hcut,
             method="silhouette",k.max=10,
             hc_func = "agnes", hc_method = "ward")

#Repeat for DIANA
hc_diana_full <- diana(state_dist_mat,diss = TRUE) %>% as.hclust()
h_d <- hc_diana_full$height
t <- 1.25
cut_d <- mean(h_d) + t * sd(h_d)
cutree(hc_diana_full, h = cut_d)
fviz_nbclust(x=state_rr_all_mat,diss=state_dist_mat,
             FUN=hcut,hc_func = "diana",
             method="silhouette",k.max=10)

#Based on this it looks like Mojena's rule is decent enough, will dynamically apply

#Nested tree for JSON
hclust_to_nested_tree <- function(hclust_obj, labels, clusters) {
  n <- length(labels)  # number of leaves
  
  get_node <- function(i) {
    if (i <= n) {
      # Leaf node (observation)
      label <- labels[i]
      return(list(
        name = label,
        cluster = unname(clusters[[label]]),
        height = 0
      ))
    } else {
      # Internal node
      merge_idx <- i - n
      left_id <- hclust_obj$merge[merge_idx, 1]
      right_id <- hclust_obj$merge[merge_idx, 2]
      
      children <- lapply(c(left_id, right_id), function(id) {
        if (id < 0) {
          get_node(-id)  # leaf (negative index)
        } else {
          get_node(id + n)  # internal node (1-based)
        }
      })
      
      return(list(
        name = paste0("node_", i),
        children = children,
        height = hclust_obj$height[merge_idx]
      ))
    }
  }
  
  get_node(2 * n - 1)  # root node index
}


#Use a dummy date to indicate the full set
FULL_SET_DATE <- as.IDate("1900-01-01")

#Snapshots for D3
state_rr_snap <- fread(paste0(fn_path_rr,"df_state_rr_snap.tsv"))
state_rr_snap <- state_rr_all %>% 
  mutate(year_date=FULL_SET_DATE) %>% 
  bind_rows(state_rr_snap)
snapshot_dates <- unique(state_rr_snap$year_date)

#Make input tibble with different date options and DIANA vs AGNES
tree_input <- expand.grid(d = snapshot_dates, di_tree = c(TRUE,FALSE)) %>% tibble()

snapshot_trees <- pmap_dfr(tree_input, function(d,di_tree) {
  snapshot_rr <- state_rr_snap %>% filter(year_date == d)
  snapshot_rr_mat <- snapshot_rr %>%
  select(c("x", "y", "RR")) %>%
  rename(name = x) %>%
  pivot_wider(names_from = y, values_from = RR) %>%
  arrange(name) %>%
  select(name, sort(setdiff(names(.), "name"))) %>%
  column_to_rownames("name") %>%
  as.matrix() %>%
  { #Remove states that have no sequence pairs
    diag_na <- is.na(diag(.))
    if (any(diag_na)) {
      keep <- !diag_na
      .[keep, keep]
      
    } else {
      .
    }
  } 
  labels <- rownames(snapshot_rr_mat)
  snapshot_rr_dist <- exp(-snapshot_rr_mat) %>% as.dist()
  

  print(snapshot_rr_mat)
  
  #ifelse will just return the first element in the object
  if(di_tree){
    snapshot_hc <- diana(snapshot_rr_dist,diss = TRUE)   
  } else{
    snapshot_hc <- agnes(snapshot_rr_dist, method = "ward", diss = TRUE)
  }
                        
  hc_as_hclust <- as.hclust(snapshot_hc)
  
  h <- hc_as_hclust$height
  t <- 1.25
  cut_h <- mean(h) + t * sd(h)
  snapshot_clusters <-cutree(hc_as_hclust,h=cut_h)
  print("Total number of clusters:")
  print(max(snapshot_clusters))
  if(di_tree){
    fviz_nbclust(x=snapshot_rr_mat,diss=snapshot_rr_dist,
                 FUN=hcut,hc_func = "diana",
                 method="silhouette",k.max=6)
  }else{
    fviz_nbclust(x=snapshot_rr_mat,diss=snapshot_rr_dist,
                 FUN=hcut,hc_func = "agnes",hc_method = "ward",
                 method="silhouette",k.max=6)
  }

  readline("Press [Enter] for next tree")
  
  
  cluster_list <- as.list(snapshot_clusters)
  
  tree <- hclust_to_nested_tree(hc_as_hclust, labels = labels, clusters = snapshot_clusters)
  
  dateStr <- ifelse(d == FULL_SET_DATE,"Full Time Period",as.character(d))
  treeMethod <- ifelse(di_tree,"DIANA","AGNES")
  c <- ifelse(di_tree,snapshot_hc$dc,snapshot_hc$ac)
  
  
  tibble(
    date = dateStr,
    coef = c,
    cluster_assignments = list(cluster_list),
    tree = list(tree),
    method = treeMethod,
    hc_tree = list(hc_as_hclust)
  )
})

write_json(snapshot_trees %>% select(-hc_tree), "tree_snapshots.json", pretty = TRUE, auto_unbox = TRUE)
meth_list <- c("DIANA","AGNES")
NUM_DATE_PAIRS <- length(snapshot_dates) - 2 #Exclude the dummy date for full time
for(m in meth_list){
  t_list <- snapshot_trees %>% filter(method==m)
  plots <- list()
  fn_plot_path <- paste0("figs/",scenario,"/tangle/",m,"/")
  dir.create(fn_plot_path,recursive = TRUE,showWarnings = FALSE)
  for(i in 1:NUM_DATE_PAIRS){
    print(i)
    old_tree <- (t_list %>% 
                   filter(date == as.character(snapshot_dates[i + 1])) %>% 
                   select(hc_tree))[[1]][[1]] # I hate lists
    new_tree <- (t_list %>% 
                   filter(date == as.character(snapshot_dates[i + 2])) %>% 
                   select(hc_tree))[[1]][[1]] # I hate lists
    
    #Limit to only shared 
    intersect_pair <- dendextend::intersect_trees(as.dendrogram(old_tree),as.dendrogram(new_tree))
    
    entangle_score <- tryCatch(
      {dendextend::entanglement(intersect_pair)},
      error = function(err){return(NA)}
    )
    cpcc <- dendextend::cor_cophenetic(intersect_pair)
    print(cpcc)
    
    jpeg(paste0(fn_plot_path,"subplot_",i,".jpg"),width=600,height=600,res=150)
    par(oma = c(0, 3, 0, 3))
    dendextend::tanglegram(
      intersect_pair,
      lab.cex = 1e-6,
      margin_inner = 0,
      main_left = snapshot_dates[i + 1],
      main_right = snapshot_dates[i + 2],
      common_subtrees_color_branches = TRUE,
      highlight_distinct_edges = FALSE,
      cex_main = 1.2,
    )
    mtext(paste0("Entanglement: ",round(entangle_score,3)), side = 1, line = 2, cex = 0.9)
    dev.off()
  }
  
  # Read images
  imgs <- image_read(sprintf("%ssubplot_%i.jpg",fn_plot_path,1:NUM_DATE_PAIRS))
  # Combine into 3x3 quilt
  quilt <- image_montage(imgs, tile = "3x3", geometry = "600x600+10+10")
  # Save the quilt
  image_write(quilt, paste0(fn_plot_path,"tanglegram_quilt.jpg"),format = "jpg")
}

NUM_DATE_PAIRS <- length(snapshot_dates) - 2 #Exclude the dummy date for full time
for(m in meth_list){
  t_list <- snapshot_trees %>% filter(method==m)
  plots <- list()
  fn_plot_path <- paste0("figs/",scenario,"/tangle/",m,"/")
  dir.create(fn_plot_path,recursive = TRUE,showWarnings = FALSE)
  for(i in 1:NUM_DATE_PAIRS){
    print(i)
    old_tree <- (t_list %>% 
                   filter(date == as.character(snapshot_dates[i + 1])) %>% 
                   select(hc_tree))[[1]][[1]] # I hate lists
    new_tree <- (t_list %>% 
                   filter(date == as.character(snapshot_dates[i + 2])) %>% 
                   select(hc_tree))[[1]][[1]] # I hate lists
    
    #Limit to only shared 
    intersect_pair <- dendextend::intersect_trees(as.dendrogram(old_tree),as.dendrogram(new_tree))
    
    entangle_score <- tryCatch(
      {dendextend::entanglement(intersect_pair)},
      error = function(err){return(NA)}
    )
    cpcc <- dendextend::cor_cophenetic(intersect_pair)
    print(cpcc)
    
    jpeg(paste0(fn_plot_path,"subplot_",i,".jpg"),width=600,height=600,res=150)
    par(oma = c(0, 3, 0, 3))
    dendextend::tanglegram(
      intersect_pair,
      lab.cex = 1e-6,
      margin_inner = 0,
      main_left = snapshot_dates[i + 1],
      main_right = snapshot_dates[i + 2],
      common_subtrees_color_branches = TRUE,
      highlight_distinct_edges = FALSE,
      cex_main = 1.2,
    )
    mtext(paste0("Entanglement: ",round(entangle_score,3)), side = 1, line = 2, cex = 0.9)
    dev.off()
  }
  
  # Read images
  imgs <- image_read(sprintf("%ssubplot_%i.jpg",fn_plot_path,1:NUM_DATE_PAIRS))
  # Combine into 3x3 quilt
  quilt <- image_montage(imgs, tile = "3x3", geometry = "600x600+10+10")
  # Save the quilt
  image_write(quilt, paste0(fn_plot_path,"tanglegram_quilt.jpg"),format = "jpg")
}


#Tseries trees will use all the time points to make continuous comparisons
#Snapshots for D3
state_rr_tseries <- fread(paste0(fn_path_rr,"df_state_rr_series.tsv"))
tseries_dates <- unique(state_rr_tseries$date)

#Make input tibble with different date options and DIANA vs AGNES
tree_input <- expand.grid(d = tseries_dates, di_tree = c(TRUE,FALSE)) %>% tibble()

tseries_trees <- pmap_dfr(tree_input, function(d,di_tree) {
  rr_t <- state_rr_tseries %>% filter(date == d)
  rr_t_mat <- rr_t %>%
    select(c("x", "y", "RR")) %>%
    rename(name = x) %>%
    pivot_wider(names_from = y, values_from = RR) %>%
    arrange(name) %>%
    select(name, sort(setdiff(names(.), "name"))) %>%
    column_to_rownames("name") %>%
    as.matrix() %>%
    { #Remove states that have no sequence pairs
      diag_na <- is.na(diag(.))
      if (any(diag_na)) {
        keep <- !diag_na
        .[keep, keep]
        
      } else {
        .
      }
    } 
  labels <- rownames(rr_t_mat)
  rr_t_dist <- exp(-rr_t_mat) %>% as.dist()
  
  
  print(rr_t_mat)
  
  #ifelse will just return the first element in the object
  if(di_tree){
    hc_t <- diana(rr_t_dist,diss = TRUE)   
  } else{
    hc_t <- agnes(rr_t_dist, method = "ward", diss = TRUE)
  }
  
  hc_as_hclust <- as.hclust(hc_t)
  
  h <- hc_as_hclust$height
  t <- 1.25
  cut_h <- mean(h) + t * sd(h)
  clusters_t <-cutree(hc_as_hclust,h=cut_h)
  print("Total number of clusters:")
  print(max(clusters_t))
  
  cluster_list <- as.list(clusters_t)
  treeMethod <- ifelse(di_tree,"DIANA","AGNES")
  c <- ifelse(di_tree,hc_t$dc,hc_t$ac)
  
  tibble(
    date = d,
    coef = c,
    cluster_assignments = list(cluster_list),
    method = treeMethod,
    hc_tree = list(hc_as_hclust)
  )
})

NUM_DATE_PAIRS <- length(tseries_dates) - 1
tree_input<-expand.grid(i = 1:NUM_DATE_PAIRS,m = meth_list)

tseries_corr <- pmap_dfr(tree_input, function(i,m){
  t_list <- tseries_trees %>% filter(method==m)
  old_tree <- (t_list %>% 
                 filter(date == tseries_dates[i]) %>% 
                 select(hc_tree))[[1]][[1]] # I hate lists
  new_tree <- (t_list %>% 
                 filter(date == tseries_dates[i+1]) %>% 
                 select(hc_tree))[[1]][[1]] # I hate lists
  #Limit to only shared 
  intersect_pair <- dendextend::intersect_trees(as.dendrogram(old_tree),as.dendrogram(new_tree))
  
  ent <- dendextend::entanglement(intersect_pair)
  cpcc <- dendextend::cor_cophenetic(intersect_pair)
  
  tibble(
    d = tseries_dates[i],
    m = m,
    cpc_cor = cpcc,
    entangle_score = ent
  )
})

#Correlation plot
cpcc_plot <- ggplot(tseries_corr,aes(x=d,y=cpc_cor,color=m)) + 
  geom_line(alpha = 0.5) + 
  geom_point(alpha = 0.5) +
  geom_smooth(lty=2) +
  labs(
    x = "Date",
    y = "Cophenetic Correlation",
    color = "Tree Method"
  ) +
  scale_x_date(date_breaks="3 months") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust=1),
        panel.spacing.x = unit(0.5,"inch")) +
  facet_wrap(vars(m),nrow=1)

ggsave(filename=paste0("figs/",scenario,"/tangle/cpcc_plot.jpg"),
       cpcc_plot,
       dpi = 150,
       units = "in",
       width = 10,
       height = 4)

#Entanglement plot
entangle_plot <- ggplot(tseries_corr,aes(x=d,y=entangle_score,color=m)) + 
  geom_line(alpha = 0.5) + 
  geom_point(alpha = 0.5) +
  geom_smooth(lty=2) +
  labs(
    x = "Date",
    y = "Entanglement Index",
    color = "Tree Method"
  ) +
  scale_x_date(date_breaks="3 months") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust=1),
        panel.spacing.x = unit(0.5,"inch")) +
  facet_wrap(vars(m),nrow=1)

ggsave(filename=paste0("figs/",scenario,"/tangle/entangle_plot.jpg"),
       entangle_plot,
       dpi = 150,
       units = "in",
       width = 10,
       height = 4)
break


set.seed(17)
#state_NMDS <- metaMDS(state_dist_mat, k=2)
state_PCoA <- cmdscale(state_dist_mat, eig = TRUE)
state_PCoA$perc_explained <- state_PCoA$eig/sum(state_PCoA$eig) * 100 #Add a percent explained
state_tSNE <- tsne(state_dist_mat,k=2) %>% data.frame() %>% magrittr::set_rownames(rownames(state_rr_all_mat))

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


state_kclust_pcoa <- kmeans(state_PCoA$points,centers = 4)$cluster %>%
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

break
### Testing code to determine which methods to use
#Agglomerative tests
METHOD_LIST <- c("single","complete","ward","average")
agnes_df <- map_dfr(METHOD_LIST, function(m) {
  # Perform agnes clustering
  state_hc <- agnes(state_dist_mat, method = m, diss = TRUE)
  # Cut tree into k clusters
  cluster_assignments <- cutree(as.hclust(state_hc), k = NUM_CLUSTS)
  state_hclust <- tibble(
    state = names(cluster_assignments),
    cluster = factor(cluster_assignments)
  )
  
  # Dendrogram prep
  state_dendro <- dendro_data(as.hclust(state_hc))
  state_dendro_stems <- state_dendro$segments
  state_dendro_leaves <- state_dendro$labels %>%
    rename(state = label) %>%
    inner_join(state_hclust, by = "state")
  subtitle_text <- paste("Agglomerative", m, "-", dist_full_name, "Distance",
                         "AC:",round(state_hc$ac,3)) 
  
  # Dendrogram plot
  fig_state_dendro <- ggplot() +
    theme_bw() +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          panel.grid = element_blank()) +
    geom_segment(data = state_dendro_stems,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_point(data = state_dendro_leaves,
               aes(x = x, y = y, color = cluster)) +
    geom_text(data = state_dendro_leaves,
              aes(x = x, y = y, label = state, color = cluster),
              hjust = -0.1, show.legend = FALSE) +
    labs(title = "State Dendrogram",
         subtitle = subtitle_text,
         color = "Cluster Number",
         y = "Cluster Distance") +
    scale_y_reverse() +
    coord_flip() +
    ylim(max(state_dendro_stems$y), -0.25 * max(state_dendro_stems$y)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  # Map plot
  map_state_hclust <- plot_usmap(regions = "states", data = state_hclust, values = "cluster") +
    labs(title = "State Map of Hierarchical Clusters",
         subtitle = subtitle_text,
         color = "Cluster Number") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom")
  
  # Combine and save
  fig_state_combo_hclust <- ggarrange(fig_state_dendro, map_state_hclust, nrow = 2)
  fn_state_combo_hclust <- paste0("figs/", scenario, "/clust/state_hclust_", dist_type, "_agg_", m, ".jpg")
  ggsave(fn_state_combo_hclust, fig_state_combo_hclust, width = 2400, height = 4800, units = "px")
  
  # Return results as a row
  tibble(
    method = m,
    agglomerative_coefficient = state_hc$ac,
    agnes_model = list(state_hc),
    cluster_assignments = list(state_hclust)
  )
})

#Divisive as an alternative
# Perform DIANA clustering
state_hc <- diana(state_dist_mat, diss = TRUE)
# Cut tree into k clusters
cluster_assignments <- cutree(as.hclust(state_hc), k = NUM_CLUSTS)
state_hclust <- tibble(
  state = names(cluster_assignments),
  cluster = factor(cluster_assignments)
)

# Dendrogram prep
state_dendro <- dendro_data(as.hclust(state_hc))
state_dendro_stems <- state_dendro$segments
state_dendro_leaves <- state_dendro$labels %>%
  rename(state = label) %>%
  inner_join(state_hclust, by = "state")
subtitle_text <- paste("Divisive Analysis", "-", dist_full_name, "Distance",
                       "DC:",round(state_hc$dc,3)) 

# Dendrogram plot
fig_state_dendro <- ggplot() +
  theme_bw() +
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank()) +
  geom_segment(data = state_dendro_stems,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = state_dendro_leaves,
             aes(x = x, y = y, color = cluster)) +
  geom_text(data = state_dendro_leaves,
            aes(x = x, y = y, label = state, color = cluster),
            hjust = -0.1, show.legend = FALSE) +
  labs(title = "State Dendrogram",
       subtitle = subtitle_text,
       color = "Cluster Number",
       y = "Cluster Distance") +
  scale_y_reverse() +
  coord_flip() +
  ylim(max(state_dendro_stems$y), -0.25 * max(state_dendro_stems$y)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")

# Map plot
map_state_hclust <- plot_usmap(regions = "states", data = state_hclust, values = "cluster") +
  labs(title = "State Map of Hierarchical Clusters",
       subtitle = subtitle_text,
       color = "Cluster Number") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")

# Combine and save
fig_state_combo_hclust <- ggarrange(fig_state_dendro, map_state_hclust, nrow = 2)
fn_state_combo_hclust <- paste0("figs/", scenario, "/clust/state_hclust_", dist_type, "_diana.jpg")
ggsave(fn_state_combo_hclust, fig_state_combo_hclust, width = 2400, height = 4800, units = "px")
