library(tidyverse)
library(data.table)
library(dbplyr)
library(duckdb)
library(ggplot2)
library(quantreg)
library(splines)
library(MASS)
library(mgcv)

select <- dplyr::select #Default select to dplyr
SCENARIO <-"USA_mini"

curr_theme <- theme_set(theme_bw())
curr_theme <- theme_update(axis.title = element_text(size=16),
                           axis.text = element_text(size=10))

fn_db <- paste0("db_files/db_",SCENARIO,".duckdb")
con <- DBI::dbConnect(duckdb(),fn_db)
df_meta <- tbl(con,"metadata")
df_pairs <- tbl(con,"pairs")

get_strain_date <- function(s){
  d <- filter(df_meta, strain == s) |> 
    select(date) |> 
    pull()
  return(as.Date(d))
}

df_pairs <- df_pairs |>
  # Join for strain_1 and clean date
  left_join(
    df_meta |> 
      filter(str_detect(date, "^\\d{4}-\\d{2}-\\d{2}$")) |> 
      mutate(date_1 = as.Date(date)) |> 
      select(strain, Nextstrain_clade, date_1), #Also import nextstrain clade along the way
    by = c("strain_1" = "strain")
  ) |>
  # Join for strain_2 and clean date
  left_join(
    df_meta |> 
      filter(str_detect(date, "^\\d{4}-\\d{2}-\\d{2}$")) |> 
      mutate(date_2 = as.Date(date)) |> 
      select(strain, date_2),
    by = c("strain_2" = "strain")
  ) |>
  # Now compute derived columns
  mutate(
    transmit_date = pmax(date_1, date_2, na.rm = TRUE),
    early_date = pmin(date_1, date_2, na.rm=TRUE),
    date_diff = abs(date_1 - date_2)
)
delta_t <- select(df_pairs,c(date_diff,transmit_date,early_date,Nextstrain_clade)) |> collect() 
mean_delta <- mean(delta_t$date_diff,na.rm=TRUE)
sd_delta <- sd(delta_t$date_diff,na.rm=TRUE)
UB_delta <- 200 #To remove extreme outliers

seq_counts <- df_meta |>
  group_by(date) |>
  summarize(count = n()) |>
  collect() |>
  mutate(date=as.Date(date))
              


hist_max <- 10 #Can look at automating adjusting this more
#Histogram of date_diff values
plot_date_diff <- ggplot(delta_t, aes(x = date_diff)) +
  # Histogram: shown on primary y-axis (left)
  geom_histogram(
    aes(y = 100 * after_stat(count) / sum(after_stat(count))),
    binwidth = 1,
    fill = "steelblue1",
    color = "black",
    boundary = 0
  ) +
  # ECDF: scaled to visually match histogram height
  stat_ecdf(
    aes(y = ..y.. * hist_max),
    geom = "line",
    color = "red",
    size = 1,
    lty =2
  ) +
  scale_x_continuous(
    name = "Difference between collection (days)",
    breaks = seq(0, 100, by = 10),
    limit = c(0,100)
  ) +
  scale_y_continuous(
    name = "Percent of pairs",  # Histogram axis
    sec.axis = sec_axis(
      trans = ~ . * (100 / hist_max),  # Rescale back to 0–100%
      name = "Cumulative %",
      breaks=seq(0,100,by=20)
    )
  ) 
ggsave(filename = paste0("figs/",SCENARIO,"/time/ddiff_hist.jpg"),
       plot_date_diff,
       width = 8,
       height = 4,
       units = "in",
       dpi = 300
)

#Delta-time over time, will put on hold for now 
#Must use numeric for splines and factor for Nextstrain_clade
# dT_gam <- gam(date_diff ~ s(as.numeric(early_date)),
#               data=delta_t)
# dT_gam_clade <- gam(date_diff ~ s(as.numeric(early_date)) + 
#                 as.factor(Nextstrain_clade),
#               data=delta_t)
# 
# #Ratio between two y-axes
# seq_delta_ratio <- max(seq_counts$count)/UB_delta
# 
# #Resolution for density calculation
# range_x <- as.numeric(range(delta_t$transmit_date))
# range_y <- range(delta_t$date_diff)
# x_days <- diff(range_x)      # total number of days on x
# y_days <- diff(range_y)      # total number of days on y
# n_x <- ceiling(x_days / 7)
# n_y <- ceiling(y_days / 1)
# 
# plot_deltatime <-  ggplot(delta_t,aes(x=transmit_date,y=date_diff)) +
#   stat_density_2d(
#     aes(fill = after_stat(log10(density + 1e-6))),
#     geom = "raster",
#     contour = FALSE,
#     n = c(n_x, n_y)
#   ) +
#   scale_fill_gradient(
#     low = "white",
#     high = "blue",
#     name = "log10(Density)"
#   ) +
#   geom_smooth(
#     method = "gam",
#     color = "red",
#     se = FALSE
#   )+
#   #geom_smooth(
#   #  method = "rlm",
#   #  formula = y ~ splines::ns(x, df = 15),
#   #  color="red",
#   #  se = FALSE
#   #)+
#   geom_line(data=seq_counts,aes(x=date,y=count/seq_delta_ratio),
#             alpha=0.4) +
#   scale_x_date(
#     name = "Date of later sequenece"
#   ) +
#   scale_y_continuous(
#     name = paste0( "\u0394","Time (days)"),
#     limit = c(0,UB_delta),
#     sec.axis = sec_axis(
#       trans = ~ . * (seq_delta_ratio),
#       name = "Sequence Counts",
#     )
#   )
# 
# 
# ggsave(filename = paste0("figs/",SCENARIO,"/time/deltat_time_plot.jpg"),
#        plot_deltatime,
#        width = 8,
#        height = 4,
#        units = "in",
#        dpi = 300
# )


#Similar to dbWriteTable, except just executes the SQL command in DuckDB
#This way I don't load the full table into memory
df_pairs %>% compute(name = "pairs_time", 
                     temporary = FALSE)

DBI::dbDisconnect(con)
