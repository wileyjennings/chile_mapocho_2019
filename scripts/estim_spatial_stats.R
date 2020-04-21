# Doc Info ----------------------------------------------------------------

# Title: Estimate statistics for Chile Mapocho samples across space
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 31 Dec 2019

# Description: Estimate summary statistics, correlations, and group comparisons
# across sampling locations.


# Dependencies ------------------------------------------------------------

required_packages <- c("astsa", "broom", "car", "dplyr", "geepack", "here", 
                       "lubridate", "olsrr", "readr")
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

# Number of field values censored for each target.
summary_cens <- 
  water %>%
  filter(campaign %in% c("spatial", "temporal", "hourly")) %>%
  group_by(target) %>%
  count(cens) %>%
  mutate(total = sum(n),
         frac_cens = n/total)

# Median, min, max and noro positive at each station.
summary_stats <- 
  water %>%
  filter(campaign %in% c("spatial", "temporal")) %>%
  group_by(target, location) %>%
  summarize(med = median(l10_100ml_cens, na.rm=T),
            min_conc = min(l10_100ml_cens),
            max_conc = max(l10_100ml_cens),
            noro_pos = sum(!near(l10_100ml_cens, log10(60))),
            n = n())


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
  unnest(tidied)

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
plot_gee_resids <- recordPlot(plot(gee_crass_hf$residuals, type = "b"))
plot_gee_autocor <- recordPlot(acf2(gee_crass_hf$residuals))  # No autocorrelation issues! 
gee_coefs <- summary(gee_crass_hf)$coefficients
# Summarize correlation.
crass_hf_cor = tibble(
  beta = gee_coefs[2, 1],
  ci_lo = gee_coefs[2, 1] - 1.96*gee_coefs[2, 2],
  ci_hi = gee_coefs[2, 1] + 1.96*gee_coefs[2, 2],
  p_val = gee_coefs[2, 4])


# Write results -----------------------------------------------------------

write_result(summary_cens)
write_result(summary_stats)
location_levene <- 
  unnest(location_clustering, tidied_levene) %>% select(-tidied_kruskal)
location_kruskal <- 
  unnest(location_clustering, tidied_kruskal) %>% select(-tidied_levene)
write_result(location_levene)
write_result(location_kruskal)
write_result(crass_hf_cor)

write_tif(plot_gee_resids, "gee_resids.tif")
write_tif(plot_gee_autocor, "gee_autocor.tif")

