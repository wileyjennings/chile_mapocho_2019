# Doc Info ----------------------------------------------------------------

# Title: Plot Chile Mapocho samples
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 27 Dec 2019

# Description: Plot spatial campaign data.


# Requirements ------------------------------------------------------------

required_packages <- c("here", "dplyr", "forcats", "ggbeeswarm", "ggplot2")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data
water <- readRDS(here::here("data", "processed", "water.rds"))


# Global settings ---------------------------------------------------------

theme_set(theme_bw())
scale_colour_discrete <- scale_colour_viridis_d


# Plot spatial campaign data ----------------------------------------------

# Filter spatial data and recode factors for plotting.
water_spatial <- 
  water %>%
  filter(target != "tc" & campaign %in% c("spatial")) %>%
  mutate(target = fct_relevel(target, "crass", "hf183", "noro", "ec", "ent"),
         target = fct_recode(target, crAssphage = "crass", HF183 = "hf183", 
                             Norovirus = "noro", `E. coli` = "ec", 
                             Enterococci = "ent"))

# Plot microbial data facetted by organism.
plot_spatial_facet_target <-
  water_spatial %>%
  filter(!is.na(detect)) %>%
  ggplot(., aes(x = location, y = l10_100ml_cens)) +
  # geom_hline(data = water_spatial %>% filter(vol_mce == 100),
  #            aes(yintercept = l10_100ml_lod), color = "gray50",
  #            linetype = 2) +
  # geom_hline(data = water_spatial %>% filter(target %in% c("ec", "ent")),
  #            aes(yintercept = l10_100ml_lod), color = "gray50",
  #            linetype = 2) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(aes(shape = detect, color = detect), cex = 2) +
  scale_shape_manual(values = c(19, 4)) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~target) +
  labs(x = "Location", y = expression("Concentration (log"[10]*"MPN or cp/100 ml)"), 
       color = "Organism") +
  theme(legend.position = "none")


# Write figures to file ---------------------------------------------------

write_tif(plot_spatial_facet_target, "spatial_facet_target.tif", "horiz")
