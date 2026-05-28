#File: age_RR_deviance.R
#Author(s): Amin Bemanian
#Description: Calculate RR deviance by age as a proxy for how structured
# versus random the mixing pattern is for each age group.
# Deviance = RMS of log10(RR) across all interaction partners.
# Under random mixing, RR=1 for all pairs so deviance=0.
# Higher deviance indicates more structured (assortative) mixing.

library(tidyverse)
library(argparse)
library(viridis)

source("scripts/color_schemes.R")

collect_args <- function(){
  parser <- ArgumentParser()
  parser$add_argument('--scenario', default = "CAM_1000", type = 'character',
                      help = 'Which scenario to perform the analysis on')
  return(parser$parse_args())
}

args <- collect_args()
scenario <- args$scenario

# RR deviance: root mean square of log10(RR)
# Measures deviation from the random mixing null (RR=1, log10(RR)=0)
rr_dev <- function(RR_list){
  log_rr <- log10(RR_list)
  d <- sqrt(sum(log_rr^2)/length(log_rr))
  return(d)
}

# --- Overall age deviance (all pairs, no geographic stratification) ---
AGE_CAP <- 80

fn_age_RR <- paste0("results/", scenario, "/df_RR_by_age_class.tsv")
df_age_RR <- read_tsv(fn_age_RR, show_col_types = FALSE) %>%
  mutate(x_age = substr(x, 1, 2) %>% as.numeric(),
         y_age = substr(y, 1, 2) %>% as.numeric()) %>%
  filter(x_age < AGE_CAP & y_age < AGE_CAP)

df_RR_dev <- df_age_RR %>%
  group_by(x) %>%
  summarize(deviance = rr_dev(RR),
            N = median(N_x),
            .groups = "drop") %>%
  mutate(age = substr(x, 1, 2) %>% as.numeric())

# Deviance by age curve
plot_age_dev <- ggplot(df_RR_dev, aes(x = age, y = deviance)) +
  geom_point() +
  geom_smooth(method = "gam") +
  theme_bw() +
  labs(x = "Age", y = "RR Deviance")

# Sensitivity: deviance vs sample size, colored by age
plot_dev_sens <- ggplot(df_RR_dev, aes(x = N, y = deviance, color = age)) +
  geom_point() +
  theme_bw() +
  scale_colour_viridis() +
  labs(x = "Number of Sequences", y = "RR Deviance", color = "Age")

fn_age_dev_plot <- paste0("figs/", scenario, "/age_RR_deviance.jpg")
fn_dev_sens_plot <- paste0("figs/", scenario, "/age_RR_deviance_sensitivity.jpg")

ggsave(plot_age_dev, filename = fn_age_dev_plot, dpi = 300, units = "in", width = 7, height = 3)
ggsave(plot_dev_sens, filename = fn_dev_sens_plot, dpi = 300, units = "in", width = 7, height = 5)
message(paste0("Overall deviance plots saved to ", fn_age_dev_plot, " and ", fn_dev_sens_plot))

# --- Geographic stratification: same state, same region, inter-region ---
fn_age_state_RR <- paste0("results/", scenario, "/df_RR_by_age_state.tsv")
df_age_state_RR <- read_tsv(fn_age_state_RR, show_col_types = FALSE) %>%
  mutate(x_age = substr(x, 1, 2) %>% as.numeric(),
         y_age = substr(y, 1, 2) %>% as.numeric()) %>%
  filter(x_age < AGE_CAP & y_age < AGE_CAP)

# Classify geographic relationship
df_age_state_RR <- df_age_state_RR %>%
  mutate(
    geo_class = case_when(
      sameState == TRUE  ~ "Same State",
      sameRegion == TRUE ~ "Same Region",
      TRUE               ~ "Inter-Region"
    )
  )

# Calculate deviance for each age × geographic class
df_geo_dev <- df_age_state_RR %>%
  group_by(x, geo_class) %>%
  summarize(deviance = rr_dev(RR),
            N = median(N_x),
            .groups = "drop") %>%
  mutate(
    age = substr(x, 1, 2) %>% as.numeric(),
    geo_class = factor(geo_class, levels = c("Same State", "Same Region", "Inter-Region"))
  )

# Combined deviance curve by geographic class
plot_geo_dev <- ggplot(df_geo_dev, aes(x = age, y = deviance, color = geo_class)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "gam") +
  theme_bw() +
  scale_color_brewer(palette = "Set1", name = "Geographic Relationship") +
  labs(x = "Age", y = "RR Deviance") +
  theme(legend.position = "bottom")

fn_geo_dev_plot <- paste0("figs/", scenario, "/age_RR_deviance_geographic.jpg")
ggsave(plot_geo_dev, filename = fn_geo_dev_plot, dpi = 300, units = "in", width = 7, height = 2.5)
ggsave(plot_geo_dev, filename = sub("\\.jpg$", ".svg", fn_geo_dev_plot), units = "in", width = 7, height = 2.5)

# Compact version sized to slot beneath the subset-heatmap row in the stitched figure
plot_geo_dev_compact <- plot_geo_dev +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))
ggsave(plot_geo_dev_compact,
       filename = paste0("figs/", scenario, "/age_RR_deviance_geographic_compact.png"),
       dpi = 300, units = "in", width = 5, height = 3.33)
ggsave(plot_geo_dev_compact,
       filename = paste0("figs/", scenario, "/age_RR_deviance_geographic_compact.svg"),
       units = "in", width = 5, height = 3.33)
message(paste0("Geographic deviance plot saved to ", fn_geo_dev_plot))

cat("\nAll deviance plots complete.\n")
