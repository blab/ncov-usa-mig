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
library(patchwork)

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
dir.create(paste0("figs/",scenario,"/age_heatmaps"), showWarnings = FALSE) #Subdir for all age heatmap outputs
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
    # Get unique age values and filter to every 5 years
    age_breaks <- unique(c(filtered_data$x, filtered_data$y)) %>%
      sort() %>%
      grep("[05]y$", ., value = TRUE)
    print(age_breaks)
    p <- p +
      scale_x_discrete(name = "Age Group", breaks = age_breaks) +
      scale_y_discrete(name = "Age Group", breaks = age_breaks) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = AXIS_SIZE * 1.25),
            axis.text.y = element_text(size = AXIS_SIZE * 1.25),
            axis.title = element_text(size = AXIS_SIZE * 1.25))
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

# --- Hardcoded 0-24yo subset (higher RR ceiling, 5-year labels) ---
# Standalone block: zooms into childhood/young-adult ages where the upper
# tail of RR is compressed in the full plot. Uses UB = 3 to expand contrast
# and labels every 5 years instead of every 10.
SUBSET_UB <- 2
SUBSET_LB <- 0.5

age_rr_young <- age_rr %>%
  filter(x <= "24y", y <= "24y") %>%
  rowwise() %>%
  mutate(fill_RR = max(min(log10(RR), log10(SUBSET_UB)), log10(SUBSET_LB))) %>%
  ungroup()

young_breaks <- unique(c(age_rr_young$x, age_rr_young$y)) %>%
  sort() %>%
  grep("[05]y$", ., value = TRUE)

# Native font size for panel B: 2x panel A's AXIS_SIZE * 1.25, offsetting the
# opposite scaling the two panels get in scripts/stitch/stitch_age.py.
YOUNG_FONT <- AXIS_SIZE * 1.25 * 2

# Highlight the 18yo row/column (and its diagonal cell) with bold black outlines.
young_levels <- sort(unique(c(age_rr_young$x, age_rr_young$y)))
HL_AGE <- "18y"
hl_idx <- match(HL_AGE, young_levels)
n_young <- length(young_levels)

age_heatmap_young <- ggplot(age_rr_young, aes(x = x, y = y, fill = fill_RR)) +
  geom_tile() +
  RR_log_grad(LB = SUBSET_LB, UB = SUBSET_UB) +
  annotate("rect",
           xmin = hl_idx - 0.5, xmax = hl_idx + 0.5,
           ymin = 0.5, ymax = n_young + 0.5,
           fill = NA, colour = "black", linewidth = 0.8) +
  annotate("rect",
           xmin = 0.5, xmax = n_young + 0.5,
           ymin = hl_idx - 0.5, ymax = hl_idx + 0.5,
           fill = NA, colour = "black", linewidth = 0.8) +
  annotate("text", x = hl_idx, y = hl_idx, label = "*",
           colour = "black", fontface = "bold", size = 10, vjust = 0.75) +
  scale_x_discrete(name = "Age Group", breaks = young_breaks) +
  scale_y_discrete(name = "Age Group", breaks = young_breaks) +
  theme_minimal() +
  # Native fonts/colourbar are set to 2x panel A's so they render at the same
  # on-page size: the stitch scales A up by 10/6 but B down by 5/6.
  guides(fill = guide_colourbar(barheight = unit(4.2, "cm"),
                                barwidth = unit(0.85, "cm"))) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = YOUNG_FONT),
        axis.text.x = element_text(angle = 45, hjust = 1, size = YOUNG_FONT),
        axis.text.y = element_text(size = YOUNG_FONT),
        legend.title = element_text(size = YOUNG_FONT * 1.1),
        legend.text = element_text(size = YOUNG_FONT * 0.9)) +
  coord_equal()

ggsave(paste0("figs/", scenario, "/age_heatmaps/young_0_24.jpg"),
       plot = age_heatmap_young,
       device = "jpeg",
       dpi = 300,
       width = 6,
       height = 6
)

ggsave(paste0("figs/", scenario, "/age_heatmaps/young_0_24.svg"),
       plot = age_heatmap_young,
       width = 6,
       height = 6,
       units = "in"
)

# Generate the three heatmaps using the function
age_heatmap_same_state <- make_age_heatmap(
  age_state_rr,
  same_state = TRUE,
  title = "Same Division",
  show_legend = FALSE
)

age_heatmap_same_region <- make_age_heatmap(
  age_state_rr,
  same_state = FALSE,
  same_region = TRUE,
  title = "Diff Div/Same Reg",
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

# Vector outputs for the stitched manuscript figure.
# Saved on a shorter (6 x 5.3) canvas: the heatmap square is width-limited by the
# right-side legend, so a full 6" tall canvas leaves vertical whitespace above and
# below. Trimming the height crops that slack. Keep height in sync with A_SRC_H_IN
# in scripts/stitch/stitch_age.py.
ggsave(paste0("figs/", scenario, "/age_heatmaps/full.svg"),
       plot = age_heatmap + theme(plot.margin = margin(1, 1, 1, 1)),
       width = 6, height = 5.3, units = "in")

# Titles are rendered as facet strips rather than plot titles so they can carry a
# filled background. Each plot gets a single-level facet_wrap(~"<label>"), which
# draws one strip bar across the top; the fill comes from GEO_CLASS_COLORS
# (scripts/color_schemes.R) so the labels match the curve colors in panel D.
labelled_subset <- function(p, label){
  p +
    labs(title = NULL) +
    facet_wrap(as.formula(paste0('~ "', label, '"'))) +
    theme(
      strip.background = element_rect(fill = GEO_CLASS_COLORS[[label]], colour = NA),
      # 11pt: the longest label ("Diff Div/Same Reg") runs to the strip edge at 12.
      strip.text = element_text(colour = "white", face = "bold", size = 11,
                                margin = margin(2, 2, 2, 2)),
      plot.margin = margin(2, 2, 2, 2)
    )
}

# Add a left spacer so the row of heatmaps aligns with the deviance plot panel
# (deviance has ~0.5" of y-axis space on its left in the stitched figure).
subsets_row <- plot_spacer() +
  labelled_subset(age_heatmap_same_state,       "Same Division") +
  labelled_subset(age_heatmap_same_region,      "Diff Div/Same Reg") +
  labelled_subset(age_heatmap_different_region, "Different Region") +
  plot_layout(ncol = 4, widths = c(0.15, 1.5, 1.5, 1.5))

ggsave(paste0("figs/", scenario, "/age_heatmaps/subsets_row.svg"),
       plot = subsets_row, width = 5, height = 2.0, units = "in")
ggsave(paste0("figs/", scenario, "/age_heatmaps/subsets_row.png"),
       plot = subsets_row, width = 5, height = 2.0, dpi = 300, units = "in")