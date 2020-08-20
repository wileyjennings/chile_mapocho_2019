# Doc Info ----------------------------------------------------------------

# Title: Estimate statistics for Chile Mapocho samples across space
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 31 Dec 2019

# Description: Estimate summary statistics, correlations, and group comparisons
# across sampling locations.


# Dependencies ------------------------------------------------------------

required_packages <- c("astsa", "broom", "car", "dplyr", "geepack", "here", 
                       "lubridate", "readr")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data
water <- readRDS(here::here("data", "processed", "water.rds"))


# Construct data frame for statistics -------------------------------------

water_wide <- 
  water %>%
  filter(campaign %in% c("hourly", "spatial", "temporal")) %>%
  pivot_wider(
    ., 
    id_cols = c(sample_id, target, l10_turb, prcp_24hr, campaign, location), 
    names_from = target, 
    values_from = l10_100ml_cens)

# Wide data frame where values are mean at each site
water_wide_mean <- 
  water %>%
  filter(campaign %in% c("spatial", "temporal")) %>%
  group_by(location, target) %>%
  summarize(location_l10mean = mean(l10_100ml_cens)) %>%
  ungroup() %>%
  pivot_wider(., id_cols = location, names_from = target, 
              values_from = location_l10mean)


# Summary stats -----------------------------------------------------------

# Number of days on which temporal sampling conducted 
temporal_sample_days <- 
  water %>% 
  filter(campaign == "temporal") %>% 
  distinct(day) %>%
  count() %>%
  pull()
print(glue("Temporal sampling performed on {temporal_sample_days} days."))

# Number of samples in which target not detected for each target.
summary_detect <- 
  water %>%
  filter(campaign %in% c("spatial", "temporal", "hourly")) %>%
  group_by(target) %>%
  count(detect) %>%
  mutate(total = sum(n),
         frac_det = n/total)

# Number of field values censored for each target.
summary_cens <- 
  water %>%
  filter(campaign %in% c("spatial", "temporal", "hourly")) %>%
  group_by(target) %>%
  count(cens) %>%
  mutate(total = sum(n),
         frac_cens = n/total)

# Median, 25th, 25th percentile across all data
summary_stats <- 
  water %>%
  filter(campaign %in% c("spatial", "temporal", "hourly")) %>%
  group_by(target) %>%
  summarize(med = quantile(l10_100ml_cens, probs = 0.5),
            `25th` = quantile(l10_100ml_cens, probs = 0.25),
            `75th` = quantile(l10_100ml_cens, probs = 0.75))

# Median, min, max and noro positive at each station.
summary_stats_loc <-
  water %>%
  filter(campaign %in% c("spatial", "temporal")) %>%
  mutate(location_group = 
           ifelse(location %in% c("S1", "S2"), "upstream", ifelse(
             location %in% c("S3", "S4"), "midstream", ifelse(
               location %in% c("S5", "S6"), "downstream", ifelse(
                 location %in% c("T1", "T2"), "temporal", NA))))) %>%
  group_by(target, location_group) %>%
  summarize(med = median(l10_100ml_cens, na.rm=T),
            min_conc = min(l10_100ml_cens),
            max_conc = max(l10_100ml_cens),
            n = n())

# Count samples in which noro was detected by location
noro_det_loc <- 
  water %>%
  filter(campaign %in% c("spatial", "temporal"),
         target == "noro") %>%
  group_by(location) %>%
  count(detect, name = "detect_noro")


# Examine spatial structure in data ---------------------------------------

# Examine clustering by site, which is important to consider for other
# statistical tests. First, examine homogeneity of variance between sampling
# location, using Levene test. This is important for interpreting results
# of group difference tests, such as Kruskal-Wallis. Then use non-parametric 
# Kruskal Wallis test to test difference by location.
location_clustering <- 
  water %>% 
  filter(campaign %in% c("spatial", "temporal")) %>%
  nest(data = -target) %>%
  mutate(
    levene_loc = map(
      data, ~leveneTest(l10_100ml_cens ~ factor(location), data = .x)),
    tidied_levene = map(levene_loc, tidy),
    kruskal_loc = map(
      data, ~kruskal.test(l10_100ml_cens ~ factor(location), data = .x)),
    tidied_kruskal = map(kruskal_loc, tidy)) %>%
  select(-c(levene_loc, kruskal_loc, data))


# Pairwise correlations ---------------------------------------------------

# Correlations between noro and other indicators across sites. 
# Since there is strong stratification by site, use site-averaged values. 
# Hourly excluded because no molecular measurements were made. 
# Must pivot the data such that noro can be compared to other targets.
cor_noro <- 
  water_wide_mean %>%
  pivot_longer(., cols = c(ec, ent, hf183, crass), names_to = "target", 
               values_to = "target_val") %>%
  nest(data = -target) %>%
  mutate(
    cor_noro = map(
      data, ~cor.test(.x$target_val, .x$noro, method = "kendall")),
    tidied = map(cor_noro, tidy)) %>%
  unnest(tidied) %>%
  select(-c(data, cor_noro))

# Correlation between HF183 and crassphage across all field samples, accounting
# for spatial clustering (note: low % censored values allows parametric 
# evaluation).
gee_crass_hf <- 
  geeglm(hf183 ~ crass, 
         data = water_wide %>% 
           filter(campaign %in% c("spatial", "temporal") & !is.na(hf183)) %>%
           arrange(location) %>%
           mutate(location = factor(location)), 
         id = location,
         family = gaussian(link = "identity"),
         corstr = "exchangeable",
         std.err = "san.se")

# Inspect residual of GEE model
if(interactive()) {
  plot(gee_crass_hf$residuals, type = "b")
  acf2(gee_crass_hf$residuals)  # No autocorrelation issues! 
}
gee_coefs <- summary(gee_crass_hf)$coefficients
# Summarize correlation.
crass_hf_cor = tibble(
  beta = gee_coefs[2, 1],
  ci_lo = gee_coefs[2, 1] - 1.96*gee_coefs[2, 2],
  ci_hi = gee_coefs[2, 1] + 1.96*gee_coefs[2, 2],
  p_val = gee_coefs[2, 4])


# Write results -----------------------------------------------------------

write_result(summary_detect)
write_result(summary_cens)
write_result(summary_stats)
write_result(summary_stats_loc)
write_result(noro_det_loc)
location_levene <- 
  unnest(location_clustering, tidied_levene) %>% select(-tidied_kruskal)
location_kruskal <- 
  unnest(location_clustering, tidied_kruskal) %>% select(-tidied_levene)
write_result(location_levene)
write_result(location_kruskal)
write_result(cor_noro)
write_result(crass_hf_cor)
