# Doc Info ----------------------------------------------------------------

# Project: Mapocho River, Chile water quality survey
# Author: Wiley Jennings | Boehm Group | Stanford 
# Date: 2 Dec 2019

# Description: Use several approaches to estimate qPCR standard equations.

# Notes: Given cleaned qPCR data, will estimate equations using a mixed model 
# approach (one equation with random effects for each plate).


# Dependencies ------------------------------------------------------------

required_packages <- c("dplyr", "ggplot2", "here", "lme4", "readr")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data
qpcr <- readRDS(here::here("data", "processed", "qpcr.rds"))


# Estimate standard equations ---------------------------------------------

# Split out standards.
standards <- qpcr %>% filter(task == "STANDARD")

# Estimate mixed standard model (assuming random intercept).
standards_mixed <- estim_std_mix(standards)
standards_mixed_summary <- summarize_std(standards_mixed)
# Note: the condition R2 for a mixed model is different from R2 for lm. It 
# assumes that explained variance includes both fixed and random effects, 
# compared to residual variance


# Estimate LOD and LOQ ----------------------------------------------------

# LOD: the theoretical LOD from Bustin 2009: 3 cp/rxn, per poisson sub-sampling
# lod <- tibble(target = c("crass", "hf183", "noro"), lod_l10 = log10(3))

lod_crass <- lod_logit(standards, "crass", .prob = 0.95)
lod_hf183 <- lod_logit(standards, "hf183", .prob = 0.95)
lod_noro <- lod_logit(standards, "noro", .prob = 0.95)

loq_inspect_crass <- lloq_pool(standards, "crass", 
                       .eff = 
                         standards_mixed_summary %>% 
                         filter(target == "crass") %>%
                         pull(effic))
loq_inspect_hf183 <- lloq_pool(standards, "hf183", 
                       .eff = 
                         standards_mixed_summary %>% 
                         filter(target == "hf183") %>%
                         pull(effic))
loq_inspect_noro <- lloq_pool(standards, "noro", 
                       .eff = 
                         standards_mixed_summary %>% 
                         filter(target == "noro") %>%
                         pull(effic))

# Based on inspection, noro lOD is used as LOQ because coefficient of variation
# never exceeds 0.35, and LOQ cannot be less than LOD
lod_loq <- tibble(
  target = c("crass", "hf183", "noro"),
  lod_l10 = c(lod_crass, lod_hf183, lod_noro),
  loq_l10 = c(loq_inspect_crass$loq_l10, loq_inspect_hf183$loq_l10, 
          lod_noro)) 
# If LOQ < LOD, assigned LOD to LOQ
lod_loq <- 
  lod_loq %>%
  mutate(loq_l10 = ifelse(loq_l10 < lod_l10, lod_l10, loq_l10))

# Add LOD and LOQ to summary of mixed standards.
standards_mixed_summary <- 
  standards_mixed_summary %>%
  left_join(., lod_loq, by = "target")

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
write_csv(standards_mixed_summary, here::here("results", "standards_summary.csv"))
write_tif(standards_plot, "standards.tif")
write_tif(loq_inspect_crass$loq_plot, "loq_crass.tif")
write_tif(loq_inspect_hf183$loq_plot, "loq_hf183.tif")
write_tif(loq_inspect_noro$loq_plot, "loq_noro.tif")
