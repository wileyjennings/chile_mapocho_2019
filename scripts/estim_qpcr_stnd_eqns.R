# Doc Info ----------------------------------------------------------------

# Project: Mapocho River, Chile water quality survey
# Author: Wiley Jennings | Boehm Group | Stanford 
# Date: 2 Dec 2019

# Description: Use several approaches to estimate qPCR standard equations.

# Notes: Given cleaned qPCR data, will estimate equations using master
# calibration (pooling across plates), plate-specific (separate equation for 
# each plate), and mixed model approach (one equation with random effects for 
# each plate).


# Dependencies ------------------------------------------------------------

required_packages <- c("dplyr", "ggplot2", "here", "lme4", "readr")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))
standards <- readRDS(here::here("data", "processed", "standards.rds"))
samples <- readRDS(here::here("data", "processed", "samples.rds"))


# Estimate standard equations ---------------------------------------------

# Estimate mixed standard model (assuming random intercept).
standards_mixed <- estim_std_mix(standards)
standards_mixed_summary <- summarize_std(standards_mixed)


# Estimate LOD and LOQ ----------------------------------------------------

# LOD: the theoretical LOD from Bustin 2009: 3 cp/rxn, per poisson sub-sampling
lod <- tibble(target = c("crass", "hf183", "noro"), lod_l10 = log10(3))

# LLOQ: use logistic regression to estimate concentration at which 50% of samples
# show amplification.
loq_crass <- lloq_logit(standards, "crass", .prob = 0.5)
loq_hf183 <- lloq_logit(standards, "hf183", .prob = 0.5)
loq_noro <- lloq_logit(standards, "noro", .prob = 0.5)

# Assign noro lloq to be log10(3), because of estimation problem and estimated
# lloq < lod = 3
loq <- tibble(target = c("crass", "hf183", "noro"),
              loq_l10 = c(loq_crass, loq_hf183, log10(3)))

# Add LOD and LOQ to summary of mixed standards.
standards_mixed_summary <- 
  standards_mixed_summary %>%
  left_join(., lod, by = "target") %>%
  left_join(., loq %>% select(loq_l10, target), by = "target")

# Plot standards. This plots all across plates.
standards_plot <- 
  standards %>%
  ggplot(., aes(x = log10(copy_num), y = ct, color = target)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, formula = y ~ x) +
  theme_bw()

# Write standards summary table and plot.
saveRDS(standards_mixed_summary, 
        here::here("data", "processed", "standards_mixed_summary.rds"))
saveRDS(standards_mixed, here::here("data", "processed", "standards_mixed.rds"))
write_csv(standards_mixed_summary, here::here("figures", "standards_summary.csv"))
write_tif_wide(standards_plot, "standards.tif")
