#File: age_heatmap.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Makes a heatmap from age RR matrix
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(viridis)

source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario
scenario <- "CAM_1000"

if(scenario ==  "CAM_1000"){
  SCALE_FACTOR <- 2.5
  RR_SIZE <- 0
  AXIS_SIZE <- 8
}else{
  SCALE_FACTOR <- 1
  RR_SIZE <- 1.5
  AXIS_SIZE <- 12
}

dir.create(paste("figs/",scenario,sep="")) #Make a directory in case it doesn't already exist
fn_rr <- paste("results/",scenario,"/df_RR_by_age_class.tsv",sep="")
age_rr <- fread(fn_rr)
fn_rr <-paste("results/",scenario,"/df_RR_by_age_state.tsv",sep="")
age_state_rr <- fread(fn_rr)

#Force a fill value to be within a certain range for display purposes
#Set as bounds for RR and transform to log for display purposes
UB <- 1.5
LB <- 0.7
fill_bound <- function(x){
  log_x <- log10(x)
  bound_x <- max(min(log_x,log10(UB)),log10(LB))
  return(bound_x)
}

# Function to create age heatmaps with geographic filters
make_age_heatmap <- function(data, same_state = NULL, same_region = NULL, title = "",
                             show_legend = TRUE, show_axis = FALSE, show_rr_labels = FALSE) {
  # Start with age filter
  filtered_data <- data %>%
    filter(x <= "80yo", y <= "80yo")

  # Apply geographic filters based on non-NULL parameters (only if columns exist)
  if (!is.null(same_state) && "sameState" %in% colnames(data)) {
    filtered_data <- filtered_data %>% filter(sameState == same_state)
  }
  if (!is.null(same_region) && "sameRegion" %in% colnames(data)) {
    filtered_data <- filtered_data %>% filter(sameRegion == same_region)
  }

  # Apply fill transformation
  filtered_data <- filtered_data %>%
    rowwise() %>%
    mutate(fill_RR = fill_bound(RR))

  # Create base plot
  p <- ggplot(filtered_data, aes(x=x, y=y, fill=fill_RR)) +
    geom_tile() +
    RR_log_grad(LB = LB,UB = UB) +
    theme_minimal() +
    labs(title = title) +
    theme(plot.title = element_text(hjust=0.5)) +
    coord_equal()

  # Add RR labels if requested
  if (show_rr_labels) {
    p <- p + geom_text_repel(aes(label=round(RR, digits=2)),
                             color = "black",
                             size = RR_SIZE,
                             bg.color = "white",
                             bg.r = 0.1,
                             force = 0)
  }

  # Configure axis display
  if (show_axis) {
    # Get unique age values and filter to every 10 years
    age_breaks <- unique(c(filtered_data$x, filtered_data$y)) %>%
      sort() %>%
      grep("*0y$", ., value = TRUE)
    print(age_breaks)
    p <- p +
      scale_x_discrete(name = "Age Groups", breaks = age_breaks) +
      scale_y_discrete(name = "Age Groups", breaks = age_breaks) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = AXIS_SIZE),
            axis.text.y = element_text(size = AXIS_SIZE))
  } else {
    p <- p + theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  }

  # Optionally hide legend
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}

age_heatmap <- make_age_heatmap(
  age_rr,
  title = NULL,
  show_legend = TRUE,
  show_axis = TRUE,
  show_rr_labels = TRUE
)

fn_age_plot <- paste0("figs/",scenario,"/age_heatmaps/full",".jpg") 
  
ggsave(fn_age_plot,
       plot=age_heatmap,
       device = "jpeg",
       dpi = 300,
       width = 6,
       height = 6
)

# Generate the three heatmaps using the function
age_heatmap_same_state <- make_age_heatmap(
  age_state_rr,
  same_state = TRUE,
  title = "Same State",
  show_legend = FALSE
)

age_heatmap_same_region <- make_age_heatmap(
  age_state_rr,
  same_state = FALSE,
  same_region = TRUE,
  title = "Same Region, Different State",
  show_legend = FALSE
)

age_heatmap_different_region <- make_age_heatmap(
  age_state_rr,
  same_region = FALSE,
  title = "Different Region",
  show_legend = FALSE
)

# Save the three geographic heatmaps
HEATMAP_WIDTH <- 3
HEATMAP_HEIGHT <- 3
HEATMAP_DPI <- 300

ggsave(paste0("figs/", scenario, "/age_heatmaps/heatmap_same_state.jpg"),
       plot = age_heatmap_same_state,
       device = "jpeg",
       dpi = HEATMAP_DPI,
       width = HEATMAP_WIDTH,
       height = HEATMAP_HEIGHT,
       units = "in")

ggsave(paste0("figs/", scenario, "/age_heatmaps/heatmap_same_region.jpg"),
       plot = age_heatmap_same_region,
       device = "jpeg",
       dpi = HEATMAP_DPI,
       width = HEATMAP_WIDTH,
       height = HEATMAP_HEIGHT,
       units = "in")

ggsave(paste0("figs/", scenario, "/age_heatmaps/heatmap_different_region.jpg"),
       plot = age_heatmap_different_region,
       device = "jpeg",
       dpi = HEATMAP_DPI,
       width = HEATMAP_WIDTH,
       height = HEATMAP_HEIGHT,
       units = "in")
