#File: find_hyper_seq.R
#Author(s): Amin Bemanian
#Date: 07/07/25
#Description: Find sequences that are hyper-represented in the pairs file
#Arguments: 
#--scenario: Scenario corresponding to data files, typically will be a geographic division (e.g. USA or Washington)

library(argparse)
library(dplyr)
library(data.table)
library(ggplot2)
library(duckdb)
library(dbplyr)
library(scales)
library(purrr)
library(patchwork)  # Or use cowplot::plot_grid()


collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', type = 'character', help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}


args <- collect_args()
scenario <- "USA"

fn_db <- paste0("db_files/db_",scenario,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)

df_meta <- tbl(con,"metadata")
df_pairs <- tbl(con,"pairs") %>%
  left_join(select(df_meta, strain, Nextstrain_clade,division),
            by =join_by(strain_1 == strain)
            )

strain_counts <- df_pairs %>% 
  group_by(strain_1) %>% 
  summarize(
    n = n(),
    Nextstrain_clade = first(Nextstrain_clade),
    division = first(division)
  ) %>% 
  arrange(desc(n)) %>%
  collect()

hist_max <- 6E5
fn_strain_plot <- paste0("figs/",scenario,"/all_strain_pair_freq.jpg")
strain_hist <- ggplot(data = strain_counts,aes(x = n)) +
  geom_histogram(fill="lightblue",color="black",bins=20) +
  labs(x = "Sequence frequency in pair set (log-scaled)",
       y = "Number of sequences") +
  stat_ecdf(
    aes(y = ..y.. * hist_max),
    geom = "line",
    color = "red",
    size = 1,
    lty =2
  ) +
  geom_line(y = 0.90 * hist_max, color = "black", lty = 2) +
  scale_x_log10(labels=label_comma()) +
  scale_y_continuous(
    labels=label_comma(),
    sec.axis = sec_axis(
      trans = ~ . * (100 / hist_max),  # Rescale back to 0–100%
      name = "Cumulative %",
      breaks=seq(0,100,by=20)
      )
    ) +
  theme_bw() 
ggsave(strain_hist,filename = fn_strain_plot,
       height = 3,
       width = 5,
       dpi = 300)


strain_by_clade <- strain_counts %>%
  group_split(Nextstrain_clade)

plot_clade_hist_ecdf <- function(df) {
  clade_name <- unique(df$Nextstrain_clade)
  
  # Create the histogram object first
  hist_plot <- ggplot(df, aes(x = n)) +
    geom_histogram(bins = 20)
  
  # Use ggplot_build on the object to extract histogram counts
  hist_max <- ggplot_build(hist_plot)$data[[1]]$count %>%
    max(na.rm = TRUE)
  
  # Build the final plot using the per-clade hist_max
  ggplot(df, aes(x = n)) +
    geom_histogram(fill = "lightblue", color = "black", bins = 20) +
    stat_ecdf(
      aes(y = ..y.. * hist_max),
      geom = "line",
      color = "red",
      size = 1,
      lty = 2
    ) +
    geom_hline(yintercept = 0.90 * hist_max, color = "black", lty = 2) +
    scale_x_log10(labels = label_comma()) +
    scale_y_continuous(
      labels = label_comma(),
      sec.axis = sec_axis(
        trans = ~ . * (100 / hist_max),
        name = "Cumulative %",
        breaks = seq(0, 100, by = 20)
      )
    ) +
    labs(
      title = paste(clade_name),
      x = "Sequence frequency",
      y = "Number of sequences"
    ) +
    theme_bw()
}

clade_plots <- map(strain_by_clade, plot_clade_hist_ecdf)
fn_clade_plots <- paste0("figs/",scenario,"/clade_strain_pair_freq.jpg")
wrap_plots(clade_plots, ncol = 10) |> ggsave(filename = fn_clade_plots,
                                             height = 15,
                                             width = 40,
                                             dpi = 300)

strain_by_state <- strain_counts %>%
  group_split(division)

plot_state_hist_ecdf <- function(df) {
  state_name <- unique(df$division)
  
  # Create the histogram object first
  hist_plot <- ggplot(df, aes(x = n)) +
    geom_histogram(bins = 20)
  
  # Use ggplot_build on the object to extract histogram counts
  hist_max <- ggplot_build(hist_plot)$data[[1]]$count %>%
    max(na.rm = TRUE)
  
  # Build the final plot using the per-clade hist_max
  ggplot(df, aes(x = n)) +
    geom_histogram(fill = "lightblue", color = "black", bins = 20) +
    stat_ecdf(
      aes(y = ..y.. * hist_max),
      geom = "line",
      color = "red",
      size = 1,
      lty = 2
    ) +
    geom_hline(yintercept = 0.90 * hist_max, color = "black", lty = 2) +
    scale_x_log10(labels = label_comma()) +
    scale_y_continuous(
      labels = label_comma(),
      sec.axis = sec_axis(
        trans = ~ . * (100 / hist_max),
        name = "Cumulative %",
        breaks = seq(0, 100, by = 20)
      )
    ) +
    labs(
      title = paste(state_name),
      x = "Sequence frequency",
      y = "Number of sequences"
    ) +
    theme_bw()
}

state_plots <- map(strain_by_state, plot_state_hist_ecdf)
fn_state_plots <- paste0("figs/",scenario,"/state_strain_pair_freq.jpg")
wrap_plots(state_plots, ncol = 10) |> ggsave(filename = fn_state_plots,
                                             height = 15,
                                             width = 40,
                                             dpi = 300)
