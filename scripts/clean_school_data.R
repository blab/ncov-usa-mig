#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(lme4)
library(splines)
library(lmerTest)

# Load raw school district data
df_raw <- read_csv("data/school_district_raw.csv")

# State abbreviation to name mapping
state_mapping <- c(
  "AL" = "Alabama", "AK" = "Alaska", "AZ" = "Arizona", "AR" = "Arkansas",
  "CA" = "California", "CO" = "Colorado", "CT" = "Connecticut", "DE" = "Delaware",
  "FL" = "Florida", "GA" = "Georgia", "HI" = "Hawaii", "ID" = "Idaho",
  "IL" = "Illinois", "IN" = "Indiana", "IA" = "Iowa", "KS" = "Kansas",
  "KY" = "Kentucky", "LA" = "Louisiana", "ME" = "Maine", "MD" = "Maryland",
  "MA" = "Massachusetts", "MI" = "Michigan", "MN" = "Minnesota", "MS" = "Mississippi",
  "MO" = "Missouri", "MT" = "Montana", "NE" = "Nebraska", "NV" = "Nevada",
  "NH" = "New Hampshire", "NJ" = "New Jersey", "NM" = "New Mexico", "NY" = "New York",
  "NC" = "North Carolina", "ND" = "North Dakota", "OH" = "Ohio", "OK" = "Oklahoma",
  "OR" = "Oregon", "PA" = "Pennsylvania", "RI" = "Rhode Island", "SC" = "South Carolina",
  "SD" = "South Dakota", "TN" = "Tennessee", "TX" = "Texas", "UT" = "Utah",
  "VT" = "Vermont", "VA" = "Virginia", "WA" = "Washington", "WV" = "West Virginia",
  "WI" = "Wisconsin", "WY" = "Wyoming", "DC" = "District of Columbia"
)

# Clean and aggregate data
df_clean <- df_raw %>%
  # Parse month column (format: YYYYmMM) into date
  mutate(
    year = as.integer(substr(month, 1, 4)),
    month_num = as.integer(substr(month, 6, nchar(month))),
    date = ymd(paste(year, month_num, "01", sep = "-"))
  ) %>%
  # Convert state abbreviations to full names
  mutate(state = state_mapping[StateAbbrev]) %>%
  # Group by state and month, calculate average shares
  group_by(state, date) %>%
  summarise(
    share_inperson = mean(share_inperson, na.rm = TRUE),
    share_hybrid = mean(share_hybrid, na.rm = TRUE),
    share_virtual = mean(share_virtual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Sort by state and date
  arrange(state, date)

# Save cleaned data
write_csv(df_clean, "data/state_school_share.csv")

message("Cleaned school data saved to data/state_school_share.csv")

logit_conversion <- function(p){
  e <- 1E-6
  log10((p+e)/(1-p+e))
}

### Look at within vs between state effects
df_model <- df_raw %>% 
  mutate(logit_inperson = logit_conversion(share_inperson)) %>%
  mutate(logit_hybrid = logit_conversion(share_hybrid)) %>%
  mutate(logit_virtual = logit_conversion(share_virtual)) %>%
  mutate(
    year = as.integer(substr(month, 1, 4)),
    month_num = as.integer(substr(month, 6, nchar(month))),
    date = ymd(paste(year, month_num, "01", sep = "-"))
  ) %>%
  mutate(NCESDistrictID = factor(NCESDistrictID))
hist(df_model$logit_virtual)
lmer_district <- lmer(logit_inperson ~ bs(date) +
                        (1|StateAbbrev) +
                        (1|NCESDistrictID),data=df_model)
summary(lmer_district)
ranova(lmer_district)
