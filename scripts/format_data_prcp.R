# Doc Info ----------------------------------------------------------------

# Project: Mapocho River, Chile water quality survey
# Author: Wiley Jennings | Boehm Group | Stanford 
# Date: 2 Dec 2019

# Description: Read and format precipitation data.

# Notes:
# There are 2 problems that must be solved to format the precipitation
# data for analysis.
# 1. Dates are given in two different formats, so they must be read as character, 
# then split on regex expression, format separately, and then bound together.
# 2. Prcp resolution in not constant: some windows are 6 hours, some are 12.
# We want to compute 24-hr cumulative precipitation relative to sample collection 
# time (which is recorded with 30 min resolution). To do this, compute the 
# average precipitation in every 30 min interval. From this, compute 24 hr
# cumulative precip with a rolling window.
# Precip dates  (momento) are provided in UTC.
# Prcp values (RRR6_Valor) are in units of mm.
# Precip source: https://climatologia.meteochile.gob.cl/application/historicos/datosDescarga/330019


# Dependencies ------------------------------------------------------------

required_packages <- c("dplyr", "here", "lubridate", "RcppRoll", "readr", "zoo")
# required_packages <- c("here", "lubridate",  "readxl", "tidyverse")
lapply(required_packages, library, character.only = T)


# Read and format ---------------------------------------------------------

prcp <- readr::read_csv(
  here::here("data", "raw", "Precip_Santiago_330019_2019_Agua6Horas.csv"),
  col_types = cols(
    CodigoNacional = col_character(),
    momento = col_character(),
    RRR6_Valor = col_double(),
    Traza_Valor = col_double()
  ))
prcp <- prcp %>% select(momento, RRR6_Valor)

# Solve two date format problem: split, format, recombine
prcp_split <- prcp %>% base::split(grepl("/", .$momento))  # splits on date format
prcp_split[[1]]$datet <- 
  as.POSIXct(prcp_split[[1]]$momento, format("%d-%m-%Y %H:%M:%S"),
             tz = "UTC")
prcp_split[[2]]$datet <-
  as.POSIXct(prcp_split[[2]]$momento, format("%d/%m/%y %H:%M"),
             tz = "UTC")
# Combine list elements (now formatted) into single tibble
prcp <- 
  bind_rows(prcp_split) %>%
  select(-momento) %>%
  rename(value = RRR6_Valor) %>%
  mutate(value = ifelse(is.na(value), 0.0, value)) 

# Filter prcp for time range relevant to this study. Add observation number.
prcp <- 
  prcp %>%
  arrange(datet) %>%
  filter(datet < "2019-06-10" & datet > "2019-05-15") %>%
  mutate(obs_num = row_number())

# Create tibble of dates for each half hour during sampling period.
datet_complete <- 
  tibble(datet = seq(prcp$datet[1], prcp$datet[nrow(prcp)], 1800))
datet_complete <- left_join(datet_complete, prcp, by = "datet")

# Create grouping variable for each measurement
datet_complete <- 
  datet_complete %>%
  mutate(obs_num = ceiling(zoo::na.approx(obs_num)))

# Compute mean over each interval, assign to each 30 min interval.
datet_complete <- 
  datet_complete %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  group_by(obs_num) %>%
  summarize(mean_val = mean(value, na.rm = T)) %>%
  left_join(datet_complete, ., by = "obs_num") 

# Compute rolling sum over previous 24 hours.
datet_complete$prcp_24hr <- 
  roll_sum(datet_complete$mean_val, n = 96L, fill = 0)  # don't know why
# window is 96 for 48 units of half hours, but it works.

# Change into Santiago time zone
datet_complete$datet <- 
  with_tz(datet_complete$datet, tzone = "America/Santiago")

# Write data
saveRDS(datet_complete, here::here("data", "processed", "prcp.rds"))
