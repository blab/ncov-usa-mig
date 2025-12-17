# File: age_pyramid_example.R
# Author(s): Amin Bemanian
# Date: 12/13/25
# Description: Creates an age pyramid with 5-year age bins using census population data,
#              highlighting two random age groups

library(dplyr)
library(ggplot2)

# US Census 2020 population estimates (in thousands)
# Source: US Census Bureau
age_counts <- data.frame(
  age_bin = c("0-4y", "5-9y", "10-14y", "15-19y", "20-24y", "25-29y",
              "30-34y", "35-39y", "40-44y", "45-49y", "50-54y", "55-59y",
              "60-64y", "65-69y", "70-74y", "75-79y", "80-84y", "85-89y", "90+"),
  population = c(19175, 20196, 21190, 21067, 21866, 23475,
                 22298, 21752, 20037, 20533, 20489, 21701,
                 20868, 17475, 14111, 10673, 6550, 4267, 2127)
)

# Set factor levels for proper ordering
age_counts <- age_counts %>%
  mutate(age_bin = factor(age_bin, levels = age_bin))

# Randomly select two age groups to highlight
set.seed(42)  # For reproducibility
available_bins <- age_counts$age_bin
highlighted_bins <- sample(available_bins, 2)
coral_bin <- as.character(highlighted_bins[1])
royal_bin <- as.character(highlighted_bins[2])

cat("Highlighting age groups:\n")
cat("  Coral Red: ", coral_bin, "\n")
cat("  Royal Blue: ", royal_bin, "\n")

# Create color mapping
age_counts <- age_counts %>%
  mutate(
    highlight = case_when(
      age_bin == coral_bin ~ "Coral Red",
      age_bin == royal_bin ~ "Royal Blue",
      TRUE ~ "Default"
    )
  )

# Create the age pyramid (horizontal bar chart)
age_pyramid <- ggplot(age_counts, aes(x = age_bin, y = population, fill = highlight)) +
  geom_col(width = 0.95) +
  coord_flip() +
  scale_fill_manual(
    values = c("Coral Red" = "#FF6F61", "Royal Blue" = "#2068aa", "Default" = "gray60"),
    name = "Highlighted Groups"
  ) +
  scale_y_continuous(
    name = "Population (thousands)",
    expand = expansion(mult = c(0, 0.05)),
    labels = scales::comma
  ) +
  scale_x_discrete(name = NULL) +
  theme_bw() +
  ggtitle("US Population by Age - Census 2020") +
  theme(
    axis.title.x = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 11),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Save the plot
dir.create("figs", showWarnings = FALSE, recursive = TRUE)
fn_pyramid <- "figs/age_pyramid_example.jpg"

ggsave(
  fn_pyramid,
  plot = age_pyramid,
  device = "jpeg",
  dpi = 192,
  width = 6,
  height = 8
)

cat("\nAge pyramid saved to:", fn_pyramid, "\n")
cat("Total US population:", format(sum(age_counts$population) * 1000, big.mark = ","), "\n")
