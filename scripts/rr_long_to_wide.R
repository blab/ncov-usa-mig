library(tidyverse)

in_fn <- "results/USA/df_RR_by_state.tsv"
out_fn <- "results/USA/df_RR_state_wide.csv"

rr_long <- read_tsv(in_fn) |>
    select(x,y,RR)

rr_wide <- rr_long |>
    rename(name = x) |>
    pivot_wider(names_from = y, values_from = RR) |>
    arrange(name)

print(head(rr_wide))

write_csv(rr_wide,out_fn)