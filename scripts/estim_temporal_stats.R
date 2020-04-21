# Doc Info ----------------------------------------------------------------

# Title: Analysis of temporal campaign data
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 2 April 2020

# Description: Examine relationships between environmenta params and 
# indicator measurements in temporal sampling campaign at sites T1 and T2.
# This script explicitly considers issues of spatial pseudoreplication and 
# examines serial correlation of residuals to ensure that model estimates are
# valid.

# Requirements ------------------------------------------------------------

required_packages <- c("astsa", "broom", "dplyr", "ggplot2", "here", "olsrr", 
                       "purrr", "readr", "tidyr")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data
water <- readRDS(here::here("data", "processed", "water.rds"))


# Construct data frame for statistics -------------------------------------

##############################################################################
# Note: Censoring limits -1 (e.g., LOD-1) used for stats
##############################################################################
water <- 
  water %>%
  mutate(l10_100ml_cens = case_when(
    is.na(cens) ~ NA_real_,
    cens == "blod" ~ log10(10^l10_100ml_lod - 1),
    cens == "bloq" ~ log10(10^l10_100ml_loq - 1),
    cens == "aloq" ~ log10(10^l10_100ml_cens + 1),
    cens == "ncen" ~ l10_100ml_cens))

# Data frame for serial correlation: average only sampling sites together
water_site_mean <- 
  water %>%
  filter(campaign == "temporal") %>%
  group_by(datet, target) %>%
  summarize(val_mean = mean(l10_100ml_cens, na.rm = T)) %>%
  pivot_wider(., id_cols = c(datet, target), names_from = target, 
              values_from = val_mean) %>%
  arrange(datet)


# Structure in hourly samples ---------------------------------------------

lm_hour <- 
  water %>% 
  filter(campaign == "hourly" & target %in% c("ec", "ent")) %>%
  nest(data = -target) %>%
  mutate(lmod_datet = map(data, ~lm(l10_100ml_cens ~ datet, data = .x)),
         tidied_datet = map(lmod_datet, tidy),
         lmod_turb = map(data, ~lm(l10_100ml_cens ~ l10_turb, data = .x)),
         tidied_turb = map(lmod_turb, tidy)) %>%
  select(-c(data, lmod_datet, lmod_turb))
lm_hour_datet <- unnest(lm_hour, tidied_datet) %>% select(-tidied_turb)
lm_hour_turb <- unnest(lm_hour, tidied_turb) %>% select(-tidied_datet)


# Temporal and spatial structure in temporal campaign ---------------------

# This model examines effect of location (T1 or T2, which are close together)
# and time of day (7:30 or 13:30) on indicator concentrations
lm_temp <- 
  water %>% 
  filter(campaign == "temporal" & target != "noro") %>%
  nest(data = -target) %>%
  mutate(lmod = map(data, ~lm(l10_100ml_cens ~ location + time, data = .x)),
         tidied = map(lmod, tidy)) %>%
  unnest(tidied) %>%
  select(-c(data, lmod))
# Summary: time of day significant for FIB but not MST markers.
# Proximate sites T1 and T2 not significant for any indicator. 


# Regression: enviro vars in temporal campaign ----------------------------

# Examine effect of precipitation and turbidity, being sure to examine 
# autorcorrelation in residuals.

# First, test for residual correlation between sites T1 and T2 to avoid 
# pseudoreplication. First, estimate regression including location as covariate,
# since this is the most conservative, i.e., if there is still correlation in
# resids between locations, then there's nothing more we can do to adjust for
# this, and pseudoreplication will artificiually inflate degrees of freedom.
lm_temp_all <-
  water %>%
  filter(campaign == "temporal" & target != "noro") %>%
  nest(data = -target) %>%
  mutate(
    lmod = map(data, ~lm(l10_100ml_cens ~ prcp_24hr + l10_turb + location, 
                         data = .x)),
    tidied = map(lmod, tidy),
    glanced = map(lmod, glance),
    augmented = map(lmod, augment),
    resids = map(lmod, resid))

# Check for correlation in residuals between two sites T1 and T2.
plot_cor_T1_T2  <- list(ggplot(), ggplot(), ggplot(), ggplot())
names(plot_cor_T1_T2) <- lm_temp_all$target
for(i in seq_along(plot_cor_T1_T2)) {
  plot_cor_T1_T2[[i]] <- 
    bind_cols(select(lm_temp_all$data[[i]], 
                     campaign, datet, datet_days, location),
              select(lm_temp_all$augmented[[i]], .resid)) %>%
      pivot_wider(., id_cols = datet, names_from = location, 
                  values_from = '.resid') %>%
      ggplot(., aes(x = T1, y = T2)) +
      geom_point() +
      geom_smooth(method = "lm", formula = "y ~ x") +
      labs(title = lm_temp_all$target[[i]])
}
# Clearly correlation between two sites even after including location as a 
# covariate in the model. Therefore, average values across T1 and T2 before 
# performing regression at these two sites.

# Test for autocorrelation in residuals: https://online.stat.psu.edu/stat510/lesson/3/3.2
# Summary: Means were computed across sites T1 and T2 because residuals
# were correlated at these two sites.
# Then, linear models not considering serial correlation, with location and 
# precipitation as predictors, were run. Ent and EC did not show autocorrelation 
# in residuals, but HF183 and crass did. Both HF183 and crass were improved by 
# AR(1) models, although qqPlots still showed poor fits in the tails.
# Since time is known precisely, and intervals between samples were not uniform,
# it was decided to use a linear model with time as a continuous covariate. 
# This resulted in non-serially correlated errors for all targets, although there
# is still evidence of non-constant variance for ec and ent.
lm_temp <-
  water %>%
  filter(campaign == "temporal" & target != "noro") %>%
  group_by(datet, target, prcp_24hr, datet_days) %>%
  summarize(mean_conc = mean(l10_100ml_cens),
            mean_turb = mean(l10_turb)) %>%
  ungroup() %>%
  nest(data = -target) %>%
  mutate(
    lmod_no_datet = map(
      data, ~lm(mean_conc ~ prcp_24hr + mean_turb, data = .x)),
    tidied_no_datet = map(lmod_no_datet, tidy),
    glanced_no_datet = map(lmod_no_datet, glance),
    augmented_no_datet = map(lmod_no_datet, augment),
    resids_no_datet = map(lmod_no_datet, resid),
    lmod_datet = map(
      data, ~lm(mean_conc ~ prcp_24hr + mean_turb + datet_days, data = .x)),
    tidied_datet = map(lmod_datet, tidy),
    glanced_datet = map(lmod_datet, glance),
    augmented_datet = map(lmod_datet, augment),
    resids_datet = map(lmod_datet, resid))

# Inspect variance inflation factors for these linear models.
vif_datet <- 
  map(lm_temp$lmod_datet, ols_vif_tol) %>% 
  set_names(lm_temp$target) %>%
  bind_rows(.id = "target")
vif_no_datet <- 
  map(lm_temp$lmod_no_datet, ols_vif_tol) %>% 
  set_names(lm_temp$target) %>%
  bind_rows(.id = "target")

# Check for serial correlation in model residuals.
if(interactive()) {
  acf2(ts_resids(lm_temp, "ec", resids_no_datet))
  acf2(ts_resids(lm_temp, "ent", resids_no_datet))
  acf2(ts_resids(lm_temp, "hf183", resids_datet))
  acf2(ts_resids(lm_temp, "crass", resids_datet))
}

# Examine final model coefficients.
lm_temp_no_datet <- 
  lm_temp %>% 
  select(target, tidied_no_datet) %>% 
  unnest(tidied_no_datet) %>%
  left_join(., 
            unnest(lm_temp, glanced_no_datet) %>% 
              select(target, adj.r.squared),
            by = "target")
lm_temp_datet <- 
  lm_temp %>% 
  select(target, tidied_datet) %>% 
  unnest(tidied_datet) %>%
  left_join(., 
            unnest(lm_temp, glanced_datet) %>% 
              select(target, adj.r.squared),
            by = "target")

# Turb correlated with morning/afternoon?
lm_turb_time <-
  water %>% 
  filter(campaign == "temporal" & target != "noro") %>%
  group_by(datet, target, prcp_24hr, datet_days, time) %>%
  summarize(mean_conc = mean(l10_100ml_cens),
            mean_turb = mean(l10_turb)) %>%
  ungroup() %>%
  distinct(datet, .keep_all = T) %>%
  lm(mean_turb ~ time, data = .) %>% 
  tidy(.)

# Turb correlated with precip?
lm_turb_precip <- 
  water %>% 
  filter(campaign == "temporal" & target != "noro") %>%
  group_by(datet, target, prcp_24hr, datet_days, time) %>%
  summarize(mean_conc = mean(l10_100ml_cens),
            mean_turb = mean(l10_turb)) %>%
  ungroup() %>%
  distinct(datet, .keep_all = T) %>%
  lm(mean_turb ~ prcp_24hr, data = .) %>% 
  tidy(.) 


# Write results -----------------------------------------------------------

write_result(lm_hour_datet)
write_result(lm_hour_turb)
write_result(vif_datet)
write_result(vif_no_datet)
write_result(lm_temp_no_datet)
write_result(lm_temp_datet)
write_result(lm_turb_time)
write_result(lm_turb_precip)
map2(.x = plot_cor_T1_T2, 
     .y = paste0("T1_T2_resids_", names(plot_cor_T1_T2), ".tif"), 
     .f = ~write_tif(.x, .y))

     