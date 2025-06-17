#File: state_time_clust.R
#Author(s): Amin Bemanian
#Date: 5/22/25
#Description: Make geographical clustering plots, including HClust Tree and
# PCoA embeddings across snapshots. Also saves a JSON of trees and TSV of 
# embeddings to be usedfor visualization with D3.
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

#To do: remove the kclustering output code
# have the output be integrated into the D3 notebook
# general cleanup

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
library(patchwork)
library(magick)

select <- dplyr::select

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  parser$add_argument('--dist', type = 'character', help = "Cluster Distance Function Reciprocal (R) or Neg-Expontential (E)", default = "E")
  return(parser$parse_args())
}

args <- collect_args()
scenario <- "CAM_1000" #args$scenario
dist_type <- args$dist #'E' or 'R'

CUSTOM_SCALE <- c("#faa", "#e66", "#ea5", "#fdb", "#fd0", 
                  "#fea", "#ad4", "#bfb", "#0fd", "#9df", 
                  "#8bf", "#aaf", "#cae", "#eae", "#fcf")

region_col <- "bea_reg"
label_title <- "BEA Regions"
get_region <- function(df) df[[region_col]]

fn_path_rr <- paste0("results/",scenario,"/time_state/")
df_regions <- read_csv("data/us_states_regions.csv") %>%
  mutate(hhs_reg = factor(hhs_reg,
                          levels = c("Boston (Region I)",
                                     "New York (Region II)",
                                     "Philadelphia (Region III)",
                                     "Atlanta (Region IV)",
                                     "Chicago (Region V)",
                                     "Dallas (Region VI)",
                                     "Kansas City (Region VII)",
                                     "Denver (Region VIII)",
                                     "San Francisco (Region IX)",
                                     "Seattle (Region X)",
                                     "Mexico",
                                     "Eastern Canada",
                                     "Western Canada"))) %>%
  mutate(census_div = factor(census_div,
                          levels = c("East North Central",
                                     "East South Central",
                                     "Middle Atlantic",
                                     "Mountain",
                                     "New England",
                                     "Pacific",
                                     "South Atlantic",
                                     "West North Central",
                                     "West South Central",
                                     "Mexico",
                                     "Eastern Canada",
                                     "Western Canada"))) %>%
  mutate(bea_reg = factor(bea_reg,
                             levels = c("Far West",
                                        "Great Lakes",
                                        "Mideast",
                                        "New England",
                                        "Plains",
                                        "Rocky Mountain",
                                        "Southeast",
                                        "Southwest",
                                        "Mexico",
                                        "Eastern Canada",
                                        "Western Canada"))) %>%
  mutate(mig_networks = factor(mig_networks))
state_rr_all <- fread(paste0(fn_path_rr,"df_state_rr_all.tsv"))

state_rr_all_mat <- state_rr_all %>% select(c("x","y","RR")) %>%
  rename(name = x) %>%
  pivot_wider(names_from = y, values_from = RR) %>%
  arrange(name) %>%
  select(name, sort(setdiff(names(.), "name"))) %>%
  column_to_rownames("name") %>%
  as.matrix()

if (dist_type == 'E') {
  state_dist_mat <- exp(-state_rr_all_mat) %>% as.matrix() %>% as.dist()
  dist_full_name <- "Negative Exponential"
} else if (dist_type == 'R') {
  state_dist_mat <- 1 / (state_rr_all_mat) %>% as.matrix() %>% as.dist()
  dist_full_name <- "Reciprocal"
} else {
  stop("Unsupported dist_type: must be 'E' or 'R'")
}

#Figuring out cluster height to cut at based on the full AGNES tree
hc_agnes_full <- agnes(state_dist_mat,diss = TRUE) %>% as.hclust()
h_a <- hc_agnes_full$height
t <- 1.25
cut_a <- mean(h_a) + t * sd(h_a)
cutree(hc_agnes_full, h = cut_a)

#Use silhouette plots to confirm reasonable height on the full set
fviz_nbclust(x=state_rr_all_mat,diss=state_dist_mat,FUN=hcut,
             method="silhouette",k.max=20,
             hc_func = "agnes", hc_method = "ward")

#Repeat for DIANA
hc_diana_full <- diana(state_dist_mat,diss = TRUE) %>% as.hclust()
h_d <- hc_diana_full$height
t <- 1.25
cut_d <- mean(h_d) + t * sd(h_d)
cutree(hc_diana_full, h = cut_d)
fviz_nbclust(x=state_rr_all_mat,diss=state_dist_mat,
             FUN=hcut,hc_func = "diana",
             method="silhouette",k.max=20)
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
snapshot_dates <- unique(state_rr_snap$date)

#Make input tibble with different date options and DIANA vs AGNES
tree_input <- expand.grid(d = snapshot_dates, di_tree = c(TRUE,FALSE)) %>% tibble()

snapshot_trees <- pmap_dfr(tree_input, function(d,di_tree) {
  if(!is.na(d)){
    snapshot_rr <- state_rr_snap %>% filter(date == d)
  }else{ #This occurs for the full set where yeardate had previously been created as "1900-01-01"
    snapshot_rr <- state_rr_snap %>% filter(is.na(date))
  }
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

dendro_data <- dendro_data(snapshot_trees$hc_tree[[12]],type="rectangle")
label_df <- dendro_data$labels %>%
  left_join(df_regions, join_by(label == state)) %>%
  mutate(region_label = get_region(.))

gg_tree_all <- ggplot() +
  geom_segment(data = dendro_data$segments,
               aes(x = y, y = x, xend = yend, yend = xend)) +
  geom_label(data = label_df,
             aes(x = y - 0.02, y = x, label = label, fill = region_label),
             color = "black", hjust = 0, size = 3, show.legend = FALSE) +
  geom_point(data = label_df,
             aes(x = y, y = x, color = region_label),
             size = 3) +
  scale_fill_manual(values = CUSTOM_SCALE) +
  scale_color_manual(values = CUSTOM_SCALE) +
  guides(color = guide_legend(title = label_title)) +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0.3))) +
  labs(title = "State Clustering by Region",
       x = "RR Negative-Exponential Distance",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 10, r = 30, b = 10, l = 10)
  )

ggsave(gg_tree_all,
       filename=paste0("figs/",scenario,"/clust/htree.jpg"),
       height = 1800,
       width = 1500,
       units = "px",
       dpi = 150)


set.seed(17)
state_PCoA_final <- cmdscale(state_dist_mat, eig = TRUE, k = 6)

state_kclust_pcoa <- kmeans(state_PCoA_final$points, centers = 6)$cluster %>%
  as_tibble(rownames = "state") %>%
  rename(cluster = value) %>%
  mutate(cluster = as.factor(cluster)) %>%
  left_join(as_tibble(state_PCoA_final$points) %>% mutate(state = rownames(state_PCoA_final$points)), by = join_by(state)) %>%
  left_join(df_regions, by = join_by(state)) %>%
  mutate(region = get_region(.))

get_hull <- function(data, comp1, comp2) {
  data %>%
    group_by(region) %>%
    slice(chull(!!sym(comp1), !!sym(comp2)))
}

hull_V1_V2 <- get_hull(state_kclust_pcoa, "V1", "V2")
hull_V1_V3 <- get_hull(state_kclust_pcoa, "V1", "V3")
hull_V1_V4 <- get_hull(state_kclust_pcoa, "V1", "V4")

make_pcoa_plot <- function(df, comp_x, comp_y, hull_data) {
  ggplot(df, aes(x = .data[[comp_x]], y = .data[[comp_y]], color = region)) +
    geom_point(show.legend = FALSE,alpha=1) +
    geom_polygon(data = hull_data, 
                 aes(x = .data[[comp_x]], y = .data[[comp_y]], group = region, fill = region),
                 alpha = 0.2, show.legend = FALSE) +
    scale_color_manual(values = CUSTOM_SCALE) +
    scale_fill_manual(values = CUSTOM_SCALE) +
    labs(x = paste("PCoA", comp_x), y = paste("PCoA", comp_y)) +
    theme_bw()
}

gg_pcoa_V1V2 <- make_pcoa_plot(state_kclust_pcoa, "V1", "V2", hull_V1_V2)
gg_pcoa_V1V3 <- make_pcoa_plot(state_kclust_pcoa, "V1", "V3", hull_V1_V3)
gg_pcoa_V1V4 <- make_pcoa_plot(state_kclust_pcoa, "V1", "V4", hull_V1_V4)

# Save plots
dir.create(paste0("figs/", scenario, "/clust"), recursive = TRUE, showWarnings = FALSE)
ggsave(gg_pcoa_V1V2, filename = paste0("figs/", scenario, "/clust/pcoa_V1V2.jpg"), width = 600, height = 600, units = "px", dpi = 150)
ggsave(gg_pcoa_V1V3, filename = paste0("figs/", scenario, "/clust/pcoa_V1V3.jpg"), width = 600, height = 600, units = "px", dpi = 150)
ggsave(gg_pcoa_V1V4, filename = paste0("figs/", scenario, "/clust/pcoa_V1V4.jpg"), width = 600, height = 600, units = "px", dpi = 150)

pcoa_column <- gg_pcoa_V1V2 / gg_pcoa_V1V3 / gg_pcoa_V1V4 + plot_layout(heights = c(1,1,1))

# Dynamically set column name for US map
df_regions_usmap <- df_regions %>%
  mutate(region_var = .data[[region_col]])

gg_regions_map <- plot_usmap(regions = "state", data = df_regions_usmap, values = "region_var") +
  scale_fill_manual(values = CUSTOM_SCALE) +
  guides(fill = "none")

# Compose and save the stitched layout
gg_combined <- (pcoa_column | plot_spacer() | gg_tree_all) + plot_layout(widths = c(1, 0.05, 3))
gg_combined <- gg_combined / gg_regions_map

ggsave(
  gg_combined,
  filename = paste0("figs/", scenario, "/clust/stitched_clust_", region_col, ".png"),
  width = 3500,
  height = 4000,
  units = "px",
  dpi = 192
)

stop("Finished Region Clustering")



##########################################
##########################################
### REMAINDER OF CODE IS HISTORICAL ######
### TO BE DELETED WITH FINAL RELEASE #####
##########################################
##########################################

df_snapshot_pcoa <- map_dfr(snapshot_dates, function(d) {
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
  
  snapshot_pcoa <- cmdscale(snapshot_rr_dist,eig=TRUE,k=6)
  snapshot_pcoa <- kmeans(snapshot_pcoa$points,centers = 6)$cluster %>%
    as_tibble(rownames="state") %>%
    rename(cluster = value) %>%
    mutate(cluster = as.factor(cluster)) %>%
    left_join(as_tibble(snapshot_pcoa$points) %>% mutate(state = rownames(snapshot_pcoa$points)),by=join_by(state))
  return(
    tibble(
      snapshot_date = rep(d,length(snapshot_pcoa$state)),
      state = snapshot_pcoa$state,
      V1 = snapshot_pcoa$V1,
      V2 = snapshot_pcoa$V2,
      V3 = snapshot_pcoa$V3,
      V4 = snapshot_pcoa$V4,
      V5 = snapshot_pcoa$V5,
      V6 = snapshot_pcoa$V6,
      cluster = snapshot_pcoa$cluster
    )
  )
})
stop("Finished PCoA Analysis")

#Entanglement code
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
    mtext(paste0("Entanglement: ",round(entangle_score,3)), side = 3, line = 2, cex = 0.8)
    dev.off()
  }
  
  # Read images
  imgs <- image_read(sprintf("%ssubplot_%i.jpg",fn_plot_path,1:NUM_DATE_PAIRS))
  # Combine into 3x3 quilt
  quilt <- image_montage(imgs, tile = "5x2", geometry = "800x800+10+10")
  # Save the quilt
  image_write(quilt, paste0(fn_plot_path,"tanglegram_quilt.jpg"),format = "jpg")
}

#Other clustering methods
NMDS_k_range <- 2:30
state_NMDS_list <- lapply(NMDS_k_range,function(k){metaMDS(state_dist_mat,k=k)})
extract_stress = function(n){return(n$stress)}
NMDS_stress<-sapply(state_NMDS_list,extract_stress)
NMDS_scree <- tibble(k = NMDS_k_range, stress = NMDS_stress)
ggplot(NMDS_scree,aes(x=k,y=stress)) + 
  geom_point() + 
  geom_line() +
  geom_vline(xintercept = 6, linetype = "dashed", color = "red", linewidth = 1.2) +
  theme_bw() +
  labs(x = "Components used", y = "NMDS Stress", title = "Non-metric Multidimensional Scaling Stress Plot")
state_NMDS_final <- metaMDS(state_dist_mat,k=6)
fviz_nbclust(scores(state_NMDS_final),
             kmeans,
             method='silhouette',
             k.max = 20
)

state_kclust_nmds <- kmeans(state_NMDS_final$points,centers = 6)$cluster %>%
  as_tibble(rownames="state") %>%
  rename(cluster = value) %>%
  mutate(cluster = as.factor(cluster)) %>%
  left_join(as_tibble(state_NMDS_final$points) %>% mutate(state = rownames(state_NMDS_final$points)),by=join_by(state)) %>%
  left_join(df_regions,by = join_by(state))

hull_nmds_data <- state_kclust_nmds %>%
  group_by(hhs_reg) %>%
  slice(chull(MDS1, MDS2))  # chull returns convex hull indices
ggplot(state_kclust_nmds,aes(x=MDS1,y=MDS2,fill=hhs_reg)) +
  #geom_point() +
  geom_label(aes(label=state),alpha=0.2) +
  #geom_polygon(data = hull_nmds_data, aes(x = MDS1, y = MDS2, group = hhs_reg, fill = hhs_reg), 
  #             alpha = 0.2,
  #             show.legend = FALSE) +
  labs(x="MDS Component #1",
       y="MDS Component #2",
       color = "Census Divisions",
       title = "MDS Values for States") +
  theme_bw()


ggplot(state_kclust_nmds,aes(x=MDS1,y=MDS2,color=cluster)) + 
  geom_text(aes(label=state)) +
  theme_bw()
ggplot(state_kclust_nmds,aes(x=MDS2,y=MDS3,color=cluster)) + 
  geom_text(aes(label=state)) +
  theme_bw()


pcoa_input <- state_PCoA_full$points
set.seed(17)
state_tsne <- tsne(X = pcoa_input,k=2,perplexity = 3,max_iter = 1E4,min_cost=1E-4) %>% 
  as_tibble() %>% 
  mutate(state = rownames(state_rr_all_mat))


state_tsne$cluster <- kmeans(state_tsne[,1:2],centers=5)$cluster %>% as.factor()


ggplot(state_tsne,aes(x=V1,y=V2,color=cluster)) + 
  geom_text(aes(label=state)) +
  labs(title = "tSNE Embedding of States",color = "Cluster") +
  theme_bw()

fviz_nbclust(state_tsne[,1:2] %>% as.matrix(),
             kmeans,
             method='silhouette',
             k.max = 20)


plot_usmap(
  regions = "states",
  data=state_tsne,
  values = "cluster"
)

df_snapshot_tSNE <- map_dfr(snapshot_dates, function(d) {
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
  
  snapshot_tsne <- tsne(snapshot_rr_dist,k=2,perplexity = 5,max_iter = 5E3,min_cost=1E-4) %>% 
    as_tibble() %>% 
    mutate(state = rownames(snapshot_rr_mat))
  snapshot_tsne$cluster <-kmeans(snapshot_tsne[,1:2],centers=5)$cluster %>% as.factor()
  return(
    tibble(
      snapshot_date = rep(d,length(snapshot_tsne$state)),
      state = snapshot_tsne$state,
      V1 = snapshot_tsne$V1,
      V2 = snapshot_tsne$V2,
      cluster = snapshot_tsne$cluster
    )
  )
})
print("Completed clustering")

ggplot(df_snapshot_tSNE %>% filter(snapshot_date != as.IDate("1900-01-01")),aes(x=V1,y=V2,color=cluster)) + 
  geom_text(aes(label=state)) +
  labs(title = "tSNE Embedding of States",color = "Cluster") +
  theme_bw() +
  facet_wrap(facets = vars(snapshot_date),nrow=2)

plot_usmap(data = df_snapshot_tSNE %>% filter(snapshot_date == ("2021-01-01")),
           region="states",
           values="cluster") +
  labs(title="t-SNE Snapshot 2021-01-01 to 2021-06-31",fill="Clusters")

plot_usmap(data = df_snapshot_tSNE %>% filter(snapshot_date == ("2021-07-01")),
           region="states",
           values="cluster") +
  labs(title="t-SNE Snapshot 2021-07-01 to 2021-12-31",fill="Clusters")

#Adjusted Rand Index analyses
ari_date_indices <- 2:9
ari_comparison_index <- 6 #For point-wise comparisons
ari_tSNE<-sapply(ari_date_indices,function(i){
  d1 <- snapshot_dates[i]
  d2 <- snapshot_dates[i+1]
  c1 <- df_snapshot_tSNE %>% filter(snapshot_date == d1) %>% pull(cluster)
  c2 <- df_snapshot_tSNE %>% filter(snapshot_date == d2) %>% pull(cluster)
  ari <- mclust::adjustedRandIndex(c1,c2)
  return(ari)
})
df_ari <- tibble(
  ori_date = snapshot_dates[ari_date_indices],
  ari=ari_tSNE,
  method=rep("t-SNE",length(ari_tSNE))
)

ari_pcoa<-sapply(ari_date_indices,function(i){
  d1 <- snapshot_dates[i]
  d2 <- snapshot_dates[i+1]
  c1 <- df_snapshot_pcoa %>% filter(snapshot_date == d1) %>% pull(cluster)
  c2 <- df_snapshot_pcoa %>% filter(snapshot_date == d2) %>% pull(cluster)
  ari <- mclust::adjustedRandIndex(c1,c2)
  return(ari)
})
df_ari <- df_ari %>% bind_rows(
  tibble(
    ori_date = snapshot_dates[ari_date_indices],
    ari=ari_pcoa,
    method=rep("PCoA",length(ari_pcoa))
  )
)

ari_hclust<-sapply(ari_date_indices,function(i){
  d1 <- snapshot_dates[i]
  d2 <- snapshot_dates[i+1]
  c1 <- snapshot_trees %>% 
    filter(method == "AGNES") %>% 
    filter(date==as.character(d1)) %>%
    pull(cluster_assignments) %>%
    unlist()
  c2 <- snapshot_trees %>% 
    filter(method == "AGNES") %>% 
    filter(date==as.character(d2)) %>%
    pull(cluster_assignments) %>%
    unlist()
  ari <- mclust::adjustedRandIndex(c1,c2)
  return(ari)
})

df_ari <- df_ari %>% bind_rows(
  tibble(
    ori_date = snapshot_dates[ari_date_indices],
    ari=ari_hclust,
    method=rep("H-Clust",length(ari_tSNE))
  )
)

ggplot(df_ari,aes(x=ori_date,y=ari,color=method)) +
  geom_line() +
  geom_point() +
  ylim(-1,1) +
  geom_hline(yintercept = 0) +
  labs(
    title = "Adjusted Rand Index Between Consecutive RR Snapshots",
    x = "Date of Earlier Snapshot",
    y = "Adjusted Rand Index",
    color = "Clustering Method"
  ) + 
  theme_bw()

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



#####DIST CODE MOVE TO A DIFFERENT FILE LATER

grav_dist <- fread("data/state_grav.tsv")
state_rr_time_dist <- left_join(state_rr_snap,grav_dist,by=join_by(x==state_x,y==state_y))
state_rr_time_dist <- left_join(state_rr_time_dist,nb_dist,by=join_by(x==state_x,y==state_y))




state_rr_time_dist <- state_rr_time_dist %>% left_join(
  df_travel_out %>% select(x,y,RR_trips),
  by=join_by(x,y)
)

df_cor_dist <- map_dfr(snapshot_dates, function(d) {
  snap_rr <- state_rr_time_dist %>% 
    filter(year_date == d & x != y)
  
  tibble(
    cor = c(
      cor(1 / snap_rr$nb_dist,     snap_rr$RR, method = "s", use = "complete.obs"),
      cor(1 / snap_rr$euclid_dist, snap_rr$RR, method = "s", use = "complete.obs"),
      cor(snap_rr$grav_wt,         snap_rr$RR, method = "s", use = "complete.obs"),
      cor(snap_rr$RR_trips,        snap_rr$RR, method = "s", use = "complete.obs")
    ),
    dist = c(
      "Inv. Rook Adjacency",
      "Inv. Euclidean Distance",
      "Gravity Weight",
      "Trip RR"
    ),
    date = d
  )
})

ggplot(df_cor_dist %>% filter(date != '1900-01-01'),
       aes(x=date,y=cor,color=dist)) + 
  geom_point() + 
  geom_line() +
  theme_bw() +
  labs(x = "Snapshot Time Period",
       y = "Spearman Correlation",
       color = "Distance/Weight")

ggplot(state_rr_time_dist %>% filter(year_date == '1900-01-01' & x != y),
       aes(x=nb_dist,y=RR)) +
  geom_point(alpha = 0.1) +
  scale_y_log10(limits=c(0.1,10)) +
  geom_smooth(method = "gam") +
  labs(x = "Rook Neighbor Order",
       y = "RR Identical Sequence")

ggplot(state_rr_time_dist %>% filter(year_date == '1900-01-01' & x != y),
       aes(x=euclid_dist,y=RR)) +
  geom_point(alpha = 0.1) +
  scale_y_log10(limits = c(0.1,10)) +
  scale_x_continuous(limits = c(0,1500)) +
  geom_smooth() +
  labs(x = "Euclidean Distance (km)",
       y = "RR Identical Sequence")

ggplot(state_rr_time_dist %>% filter(year_date == '1900-01-01' & x != y),
       aes(x=grav_wt,y=RR)) +
  geom_point(alpha = 0.1) +
  scale_y_log10(limits=c(0.1,10)) +
  scale_x_log10() +
  geom_smooth() +
  labs(x = "Gravity Weight",
       y = "RR Identical Sequence")

ggplot(state_rr_time_dist %>% filter(year_date == '1900-01-01' & x != y),
       aes(x=RR_trips,y=RR)) +
  geom_point(alpha = 0.1) +
  scale_y_log10(limits=c(0.1,10)) +
  scale_x_log10() +
  geom_smooth() +
  labs(x = "RR Trip",
       y = "RR Identical Sequence")

