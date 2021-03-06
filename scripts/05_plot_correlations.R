# Doc Info ----------------------------------------------------------------

# Title: Plot Chile Mapocho samples
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 27 Dec 2019

# Description: Plot correlations of microbial and env measurements.


# Requirements ------------------------------------------------------------

required_packages <- c("here", "dplyr", "forcats", "GGally", "ggplot2", "tidyr")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data
water <- readRDS(here::here("data", "processed", "water.rds"))


# Global settings ---------------------------------------------------------

theme_set(theme_bw())
scale_colour_discrete <- scale_colour_viridis_d


# Indicator correlations --------------------------------------------------

# Wide data frame required for correlations.
water_wide <- 
  water %>%
  filter(campaign %in% c("hourly", "spatial", "temporal")) %>%
  pivot_wider(., id_cols = c(sample_id, target, l10_turb, 
                             prcp_24hr, location), 
              names_from = target, values_from = l10_100ml_cens)

plot_cor_noro <- 
  water_wide %>%
  filter(!is.na(location)) %>%
  pivot_longer(., cols = c(ec:crass), names_to = "indicator", 
               values_to = "indicator_conc") %>%
  mutate(indicator = fct_relevel(indicator, c("crass", "hf183", "ec", "ent")),
         indicator = fct_recode(indicator, crAssphage = "crass", HF183 = "hf183", 
                             `E. coli` = "ec", Enterococci = "ent")) %>%
  filter(!is.na(noro)) %>%
  ggplot(., aes(x = noro, y = indicator_conc)) +
  facet_wrap(~indicator) +
  geom_smooth(method = "lm", formula = 'y ~ x', color = "gray50") +
  geom_point(aes(color = location)) +
  scale_color_manual(values = viridis::viridis_pal(
    end = 0.90, option = "C")(8)) +
  labs(x = "Norovirus concentration (log10)", 
       y = "Fecal indicator concentration (log10)")

plot_cors_matrix <- 
  water_wide %>%
  filter(complete.cases(.)) %>%
  ggpairs(., columns = c(2:3, 5:8))


# Write figures to file ---------------------------------------------------

write_tif(plot_cor_noro, "correlation_facet_noro.tif", "horiz")
write_tif(plot_cors_matrix, "correlation_matrix.tif", "horiz")
