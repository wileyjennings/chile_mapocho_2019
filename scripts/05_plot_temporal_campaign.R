# Doc Info ----------------------------------------------------------------

# Title: Plot Chile Mapocho samples
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 27 Dec 2019

# Description: Plot temporal campaign data.


# Requirements ------------------------------------------------------------

required_packages <- c("here", "dplyr", "forcats", "ggplot2", "tidyr", "viridis")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data
water <- readRDS(here::here("data", "processed", "water.rds"))


# Global settings ---------------------------------------------------------

theme_set(theme_bw())
scale_colour_discrete <- scale_colour_viridis_d


# Format data for facet plots ---------------------------------------------

# Format temporal for plotting organism + turbidity + precip in separate plots. 
# Pivot concentration data (microbial, turbidity, and precip) to long, with
# indicator column `meas_class` to distinguish organism from turbidity from
# precip `values`.
water_meas_class <- 
  water %>%
  filter(campaign %in% c("temporal", "hourly")) %>%
  pivot_longer(., cols = c(l10_100ml_cens, l10_turb, prcp_24hr),
               names_to = "meas_class", values_to = "value") %>%
  mutate(
    meas_class = forcats::fct_recode(
      meas_class, Organism = "l10_100ml_cens", Turbidity = "l10_turb",
      Precipitation = "prcp_24hr"),
    meas_class = forcats::fct_relevel(
      meas_class, "Organism", "Turbidity", "Precipitation"))

# Since concentrations at sites T1 and T2 are highly correlated, format data
# for plotting means and SEs across the two sites.
water_meas_class_mean <- 
  water_meas_class %>%
  filter(campaign == "temporal") %>%  # only makes sense for temporal campaign
  group_by(datet, time, target, meas_class) %>%
  summarize(val_mean = mean(value, na.rm = T),
            val_n = n(),
            val_sd = sd(value, na.rm = T),
            val_se = val_sd/sqrt(val_n))


# Plot temporal campaign --------------------------------------------------

# Plot microbial concentrations by day of week.
plot_temp_weekday <- 
  water %>%
  filter(target != "tc", campaign == "temporal") %>%
  ggplot(., aes(x = weekday, y = l10_100ml_cens, color = target)) +
  geom_boxplot() +
  labs(x = "Weekday", y = "Concentration, log10(unit/100 ml)", color = "Target")

# Plot microbial concentrations split by time and faceted by location.
plot_temp_facet_location <- 
  water %>%
  filter(target != "tc", campaign == "temporal") %>%
  ggplot(., aes(x = time, y = l10_100ml_cens, 
                color = target)) +
  geom_boxplot() +
  facet_wrap(~location) +
  labs(x = "Time", y = "Concentration, log10(unit/100 ml)", color = "Target")

# Plot faceted daily means of each indicator + turbidity + precipitation.
plot_temp_mean_env <- 
  water_meas_class_mean %>%
  ggplot(., aes(x = datet, y = val_mean)) +
  geom_point(data = filter(water_meas_class_mean, 
                           meas_class == "Organism",
                           target != "noro"), 
             aes(color = target)) +
  geom_errorbar(data = filter(water_meas_class_mean, 
                              meas_class == "Organism",
                              target != "noro"), 
                aes(color = target, ymin = val_mean - val_se, 
                    ymax = val_mean + val_se),
                width = 10000) +
  geom_errorbar(data = filter(water_meas_class_mean,
                              meas_class == "Turbidity"),
                aes(ymin = val_mean - val_se, ymax = val_mean + val_se),
                color = "gray25", width = 10000) +
  geom_point(data = filter(water_meas_class_mean, meas_class == "Turbidity"), 
             aes(fill = time), shape = 21) +
  geom_col(data = filter(water_meas_class_mean, meas_class == "Precipitation"),
           position = "dodge", fill = "gray25") 
# Pretty it up!
plot_temp_mean_env <- 
  plot_temp_mean_env %>% 
  facet_env(., .nrow = 3, .num_targ = 4)


# Plot hourly sampling campaign -------------------------------------------

plot_hourly_turb <- 
  water_meas_class %>%
  filter(meas_class %in% c("Organism", "Turbidity")) %>%
  ggplot(., aes(x = time, y = value)) +
  geom_point(data = filter(water_meas_class, 
                           meas_class == "Organism",
                           campaign == "hourly"), 
             aes(color = target)) +
  geom_point(data = filter(water_meas_class, 
                           meas_class == "Turbidity",
                           campaign == "hourly"))
# Pretty it up!
plot_hourly_turb <- 
  plot_hourly_turb %>%
  facet_env(., .nrow = 2, .num_targ = 2)


# Write figures to file ---------------------------------------------------

write_tif(plot_temp_weekday, "temp_weekday.tif", "horiz")
write_tif(plot_temp_facet_location, "temp_facet_location.tif", "horiz")
write_tif(plot_temp_mean_env, "temp_env_mean.tif", "vert")
write_tif(plot_hourly_turb, "hourly_turb.tif", "horiz")

