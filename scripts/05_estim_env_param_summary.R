# Doc Info ----------------------------------------------------------------

# Title: Environmental parameter summary
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 2 April 2020

# Description: 
# Summarize enviromental parameters: precipitation, turbidity, and water temp

# Requirements ------------------------------------------------------------

requiredPackages <- c( "here", "dplyr", "readr")
lapply(requiredPackages, library, character.only = T)

# Processed data
water <- readRDS(here::here("data", "processed", "water.rds"))


# Summary stats -----------------------------------------------------------

env_table_spatial <- 
  water %>%
  filter(campaign %in% c("spatial")) %>%
  distinct(sample_id, .keep_all = T) %>%
  group_by(location) %>%
  summarize_at(c("l10_turb","temp", "prcp_24hr"), 
               list(~mean(.), ~sd(.)), na.rm = T)

env_table_temp <- 
  water %>%
  filter(campaign %in% c("temporal")) %>%
  distinct(sample_id, .keep_all = T) %>%
  group_by(location, time) %>%
  summarize_at(c("l10_turb", "prcp_24hr"), list(~mean(.), ~sd(.)), na.rm = T)

env_table_temp_aggreg_time <- 
  water %>%
  filter(campaign %in% c("temporal")) %>%
  distinct(sample_id, .keep_all = T) %>%
  group_by(location) %>%
  summarize_at(c("l10_turb", "prcp_24hr"), list(~mean(.), ~sd(.)), na.rm = T)

write_csv(env_table_spatial, here::here("results", "env_var_spatial_summary.csv"))
write_csv(env_table_temp, here::here("results", "env_var_temporal_summary.csv"))
write_csv(env_table_temp_aggreg_time, here::here("results", "env_var_temporal_aggreg_time_summary.csv"))
