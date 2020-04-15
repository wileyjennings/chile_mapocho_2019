# Doc Info ----------------------------------------------------------------

# Project: Mapocho River, Chile water quality survey
# Author: Wiley Jennings | Boehm Group | Stanford 
# Date: 2 Dec 2019

# Description: Estimate qPCR reaction concentrations, LOD and LLOQ, and 
# check for inhibition.

# Notes: Given cleaned qPCR data, standard equations and standard equation
# summaries, writes .rds of estimated sample qpcr rxn concentrations and a
# .csv inhibition summary.


# Dependencies ------------------------------------------------------------

required_packages <- c("dplyr", "here", "readr", "tidyr")
lapply(required_packages, library, character.only = T)

# Processed data
qpcr <- readRDS(here::here("data", "processed","qpcr.rds"))
standards_summary <- readRDS(
  here::here("data", "processed", "standards_mixed_summary.rds"))
standards_eqns <- readRDS(here::here("data", "processed","standards_mixed.rds"))


# Estimate qpcr rxn concentrations -------------------------------------------

# Split out qPCR unknown samples.
samples <- qpcr %>% filter(task != "STANDARD")

# Join standard curves to unknowns data. Key: plate_names + target.
# Pull out random effects
standards_raneff <-
  standards_eqns %>%
  unnest(re) %>%
  select(target, plate_name = grp, intercept_dif = condval) %>%
  left_join(., standards_summary, by = "target")
samples <- 
  standards_raneff %>%
  mutate(plate_name = as.character(plate_name)) %>%
  left_join(samples, ., by = c("plate_name", "target"))

# Compute mean replicate ct values and compute SDs of unknowns. 
# Discard SDs computed by Step One Plus software.
# Discard replicates now that computed summary values using `distinct()`.
samples <- 
  samples %>%
  group_by(plate_name, target, name, diln) %>%
  summarize(ct_mean = mean(ct, na.rm = T)) %>%
  left_join(samples, ., by = c("plate_name", "target", "name", "diln")) 

# Estimate concentrations from plate-specific calibration equations.
samples <- 
  samples %>%
  mutate(l10_cn = (ct - (intercept + intercept_dif)) / slope_l10cn)

# Replace concentrations with average and sd of each qPCR replicate
samples <- 
  samples %>%
  group_by(plate_name, target, name, diln) %>%
  summarize(l10_cn_sd = sd(l10_cn, na.rm = T),
            l10_cn = mean(l10_cn, na.rm = T)) %>% 
  left_join(x = samples %>% 
              distinct(plate_name, target, name, diln, .keep_all = T) %>%
              select(-l10_cn), 
            y = .,
            by = c("plate_name", "target", "name", "diln"))


# Check for inhibition (already checked visually) -------------------------

inhibition <- 
  samples %>%
  pivot_wider(., id_cols = c(plate_name, name, target), 
              names_from = diln, values_from = ct_mean, 
              names_prefix = "diln_") %>%
  mutate(delta_ct = diln_10 - diln_1) %>%
  filter(is.nan(diln_10) | !is.na(diln_10))
inhibition_summary <- 
  inhibition %>%
  group_by(target) %>%
  count(inhibited = delta_ct < 2.32) %>%
  mutate(inhibited = ifelse(is.na(inhibited), "Did not amplify", inhibited))

# Keep only undiluted extracts for analysis now that inhibition inspected.
samples <- samples %>% filter(diln == 1)

# Add indicator columns for LOD and LOQ.
samples <- 
  samples %>%
  mutate(cens = ifelse(
    l10_cn < lod_l10 | is.na(l10_cn), "blod", ifelse(
      l10_cn < loq_l10, "bloq", "ncen")))

# Add censored concentrations, with value at censoring limit, which is useful
# for plotting though not statistics.
samples <- 
  samples %>%
  mutate(l10_cn_cens = ifelse(cens == "blod", lod_l10, ifelse(
             cens == "bloq", loq_l10, l10_cn)))

# Write results and processed data.
saveRDS(samples, here::here("data", "processed", "samples.rds"))
write_csv(inhibition_summary, here::here("results", "inhibition_summary.csv"))
