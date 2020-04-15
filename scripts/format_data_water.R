# Doc Info ----------------------------------------------------------------

# Project: Mapocho River, Chile water quality survey
# Author: Wiley Jennings | Boehm Group | Stanford 
# Date: 2 Dec 2019

# Description: Combine all data to produce `water` df for final analysis and
# compute concentrations of qpcr targets per 100 ml water.

# Notes: Given cleaned lab, precipitation, and qpcr sample data, combines
# them into one data frame and computes qpcr concentrations per 100 ml. Then
# writes combined `water` df to .rds and writes tables to summarize numbe
# of field samples collected and concentrations of markers in sewage.


# Dependencies ------------------------------------------------------------

required_packages <- c("dplyr", "here", "lubridate", "readr")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data.
lab <- readRDS(here::here("data", "processed", "lab.rds"))
prcp <- readRDS(here::here("data", "processed", "prcp.rds"))
samples <- readRDS(here::here("data", "processed", "samples.rds"))


# Link qpcr data to prcp, field, and lab data --------------------------------

# Add precipitation data to lab data.
lab <- 
  prcp %>% 
  select(datet, prcp_24hr) %>%
  left_join(lab, ., by = "datet")  

# Combine lab data and qpcr data to yield final analysis df `water`.
# Challenges:
# 1. samples$sample_id (called `name`) is a subset oflab$sample_id because 
# not all samples were tested by qpcr (e.g., the hourly samples).
# 2. Both samples and lab contain microbial data we would like to merge in long
# format with a `method` indicator column (qpcr or IDEXX).
# 3. lab contains 95% CI of MPN estimates per 100 ml water for IDEXX, 
# while samples does not contain 95% CI estimates - it contains SDs of qpcr 
# technical replicates not adjusted to 100 ml water.
# 3. lab contains blanks and other QC data we don't care about now.

# First, row bind samples and lab: select columns from samples that want to
# keep; select matching columns from lab as well as 95% CI columns which do
# not match because want to have these in final df
water <- 
  bind_rows(samples %>% 
              select(sample_id = name, target, method, lod_l10, loq_l10, 
                     cens, l10_cn_cens, l10_cn_sd),
            lab %>% 
              select(sample_id, target, method, l10_100ml_lod = lod_l10, 
                     l10_100ml_loq = loq_l10, cens, l10_100ml_cens, 
                     l10_100ml_hi95, l10_100ml_lo95))

# Second, left_join lab  to water to add other water quality and prcp data.
water <- 
  lab %>% 
  distinct(sample_id, sample_name, campaign, temp, cond, l10_turb, vol_mce, 
           notes, datet, prcp_24hr) %>%
  left_join(water, ., by = "sample_id")

# Add day of week.
water <- water %>% mutate(weekday = wday(datet, label = T)) 


# Calculate qpcr target concentrations per 100 ml water -------------------

# Important note: do not use this function to propagate sd of concentration
# estimates. Since we are working with a log-normally distributed random 
# variable, and staying on the log-scale, the sd of the log-transformed
# data does not change when adjusting concentrations for lab processing
# variables. This is because multiplication on linear scale equals addition
# on log scale, and addition of a constant does not affect sd.
water <- calc_qpcr_water(water, l10_100ml_cens, l10_cn_cens)
water <- calc_qpcr_water(water, l10_100ml_lod, lod_l10)
water <- calc_qpcr_water(water, l10_100ml_loq, loq_l10)


# Finish formatting combined df -------------------------------------------

# Split sample name into new vars indicating location, day, and time.
water <- 
  water %>% 
  filter(campaign %in% c("temporal", "hourly", "spatial")) %>% 
  separate(sample_name, sep = "-", into = c("location", "day", "time"),
           fill = "right", remove = F) %>%
  left_join(water, ., by = names(water))
water <- 
  water %>% 
  arrange(campaign, datet, location) %>% 
  mutate(location = ifelse(grepl("h", .$sample_id), "STAN", location))

# Rename sampling locations as T1, T2, or S1-S6.
locations <- c("T1", rep("T2", 11), paste0("S", c(1:6)))
names(locations) <- c("CONC", "STAN", paste0("h", c(1:10)),
                      "A", "B", "C", "D", "E", "F")
water$location <- locations[water$location]

# Create date variable as numeric in units of days
water$datet_days <- 
  as.numeric(water$datet) / (60*60*24) - 
  min(as.numeric(water$datet) / (60*60*24), na.rm = T)

# Remove total coliform and reorder levels of target
water <- water %>% filter(target != "tc")
water$target <- factor(water$target, 
                       levels = c("ec", "ent", "hf183", "crass", "noro"))

# Summarize number of field samples collected.
num_field_samp <- 
  water %>%
  filter(campaign %in% c("temporal", "hourly", "spatial")) %>%
  group_by(target) %>%
  count()  # 79 samples of each target

# Summarize sewage marker concentrations.
sewage_conc <- 
  water %>% 
  filter(sample_id == "sewage", target %in% c("hf183", "crass")) %>% 
  select(sample_id, target, l10_100ml_cens, l10_cn_sd)


# Write processed data and summary tables ---------------------------------

saveRDS(water, here::here("data", "processed", "water.rds"))
write_csv(num_field_samp, here::here("results", "num_field_samp.csv"))
write_csv(sewage_conc, here::here("results", "sewage_conc.csv"))