# Doc Info ----------------------------------------------------------------

# Title: Examine historical FIB data for Mapocho
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 4 Jan 2019

# Description: Explore trends in fecal coliform data from Aguas Andinas.


# Dependencies ------------------------------------------------------------

requiredPackages <- c("dplyr", "ggplot2", "here", "lubridate", "purrr", "readr",
                      "readxl", "tidyr", "viridisLite")
lapply(requiredPackages, library, character.only = T)
source(here::here("scripts", "util.R"))

# Processed data
water <- readRDS(here::here("data", "processed", "water.rds"))


# Global settings ---------------------------------------------------------

theme_set(theme_bw())
scale_colour_discrete <- scale_colour_viridis_d
options(ggplot2.continuous.colour = "viridis")
options(ggplot2.continuous.fill = "viridis")
options(max.print = 10000)


# Read and format Aguas Andinas data ----------------------------------------

# Data exist in two different Excel sheets, one in long and one in wide format.
# Note: language is Spanish.
aa1 <- 
  read_excel(
    here::here("data", "raw", "fecal_coliform_compiled_Aguas_Andinas.xlsx"),
    sheet = "2005-2008", col_types = c(rep("text", 4)), na = c("NA", ""))
aa2 <- 
  read_excel(
    here::here("data", "raw", "fecal_coliform_compiled_Aguas_Andinas.xlsx"),
    sheet = "2009-2018", col_types = c(rep("text", 14)), 
    na = c("NA", "", "E", "U", "-"))
aa2 <- pivot_longer(aa2, cols = Enero:Diciembre, names_to = "month", 
                    values_to = "fecal_col_MPN100ml")
aa <- bind_rows(aa1, aa2)

# Clean and format data.
# Remove `<` since will not include censored data analysis here.
aa <- map_dfr(.x = aa, ~gsub("<", "", x = .x))
# Use decimal instead of commas to indicate fractional values.
aa <- map_dfr(.x = aa, ~gsub(",", ".", x = .x))
# Rename, convert types, and relevel columns.
aa <- aa %>% rename(fc = "fecal_col_MPN100ml")
dict_month <- 
  list("Enero" = "1", "Febrero" = "2", "Marzo" = "3", "Abril" = "4", 
       "Mayo" = "5", "Junio" = "6", "Julio" = "7", "Agosto" = "8", 
       "Septiembre" = "9", "Octubre" = "10", "Noviembre" = "11", 
       "Diciembre" = "12")
aa$month = unlist(dict_month[aa$month])
aa <- 
  aa %>% 
  mutate(date = make_date(aa$year, aa$month, "15"),
         year = as.numeric(year),
         month = as.numeric(month),
         fc = as.numeric(fc),
         l10_fc = log10(fc)) 

# Consolidate station names, since there are many stations which are labeled
# with multiple names (similar but slightly different names).
dict_stations <- as.list(
  c("Pte Lo Curro", "Pte Suecia", "Pte Bulnes/Panamericana", "Pte Pudahuel",
    "Mapocho en Ruta 68", "Farfana aguas abajo", "Rinconada de Maipu", 
    "Pte Esperanza B", "Pte Pelvin A", "Pte Pelvin B", 
    "Pte Ferroviario Talagante A", "Pte Ferroviario Talagante B", 
    "Pte El Monte A", "Pte El Monte B", "Pte Esperanza A", "Pte Esperanza M",
    "Antes de la descarga PTAS Talagante A", 
    "Antes de la descarga PTAS Talagante B", "Pte La Dehesa", 
    "Pte Bulnes/Panamericana", "Pte Vespucio", "Estero Lampa", "Pte Mapocho", 
    "Farfana aguas abajo", "Rinconada de Maipu", "Pte Esperanza", "Pte Ruta 78", 
    "Pte La Dehesa", "Pte Lo Curro", "Pte Suecia", "Pte Bulnes/Panamericana", 
    "Pte Vespucio", "Estero Lampa", "Mapocho en Ruta 68", "Farfana aguas abajo", 
    "Rinconada de Maipu", "Pte Esperanza", "Pte Pelvin A", "Pte Pelvin B", 
    "Pte Ferroviario Talagante A", "Pte Ferroviario Talagante B", 
    "Mapocho en Ruta 68", "Pte Lo Curro", "Pte Suecia", "Pte Bulnes/Panamericana",
    "Pte Pudahuel", "Mapocho en Ruta 68", "Farfana aguas abajo",
    "Rinconada de Maipu", "Pte Esperanza M", "Pte Pelvin B", 
    "Pte Ferroviario Talagante B", "Pte El Monte B"))
names(dict_stations) <- aa %>% distinct(station) %>% pull(.)
aa$station <- unlist(dict_stations[aa$station])


# Plot historical data ----------------------------------------------------

plot_hist_data_station <- 
  aa %>%
  filter(station %in% c("Pte Suecia", "Pte La Dehesa", "Pte Lo Curro", 
                        "Pte Bulnes/Panamericana", "Farfana aguas abajo", 
                        "Rinconada de Maipu")) %>%
  ggplot(., aes(x = date, y = l10_fc, color = station)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x", se = F) +
  labs(title = "Historical data from AA", y = "Fecal Coliform \nlog10 MPN/100ml",
       x = "Date")

plot_hist_data <- 
  aa %>%
  filter(station %in% c("Pte Suecia", "Pte La Dehesa", "Pte Lo Curro", 
                        "Pte Bulnes/Panamericana", "Farfana aguas abajo", 
                        "Rinconada de Maipu")) %>%
  ggplot(., aes(x = date, y = l10_fc)) +
  geom_point(color = viridis::viridis(1, begin = 0.2, end = 0.2)) +
  geom_smooth(method = "lm", formula = "y ~ x", se = F,
              color = viridis::viridis(1, begin = 0.2, end = 0.2)) +
  scale_color_viridis_d(begin = 0.2, end = 0.95) +
  labs(y = expression("Fecal coliform (log"["10"]*" MPN/100 ml)"), 
       x = "Date", color = "Data source")

plot_hist_and_our_data <- 
  bind_rows(
    water %>% 
      filter(target == "ec" & !is.na(location)) %>%
      select(Station = location, date = datet, l10_val = l10_100ml_cens) %>%
      mutate(Station = "This study",
             date = as_date(date)),
    aa %>%
      select(Station = station, date, l10_val = l10_fc) %>% 
      filter(Station %in% c("Pte Suecia", "Pte La Dehesa", 
                            "Pte Bulnes/Panamericana", "Farfana aguas abajo", 
                            "Rinconada de Maipu")) %>%
      mutate(Station = "Historical")) %>%
  ggplot(., aes(x = date, y = l10_val, color = Station)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x", se = F) +
  scale_color_viridis_d(begin = 0.2, end = 0.95) +
  labs(y = expression("FIB concentration (log"["10"]*" MPN/100 ml)"), 
       x = "Date", color = "Data source")


plot_AA_monthly_trends <- 
  aa %>%
  filter(station %in% c("Pte Suecia", "Pte La Dehesa", "Pte Bulnes/Panamericana", 
                        "Farfana aguas abajo", "Rinconada de Maipu")) %>% 
  ggplot(., aes(x = month, y = l10_fc)) +
  geom_point() +
  facet_wrap(~year) +
  scale_x_continuous(breaks = seq(1, 12, 3), 
                     labels = month(seq(1, 12, 3), label = T)) +
  geom_smooth(method = 'loess', formula = 'y ~ x') +
  labs(x = "Month", y = "Fecal coliform (log10 MPN/100 ml)")

# Seasonal trends?
plot_monthly_trend_2016to18 <- 
  aa %>%
  filter(station %in% c("Pte Suecia", "Pte La Dehesa", "Pte Bulnes/Panamericana", 
                        "Farfana aguas abajo", "Rinconada de Maipu"),
         year %in% c(2016, 2017, 2018)) %>% 
  ggplot(., aes(x = as.numeric(month), y = l10_fc)) +
  geom_point() +
  geom_smooth(method = 'loess', formula = 'y ~ x') +
  scale_x_continuous(breaks = seq(1, 12, 3), 
                     labels = month(seq(1, 12, 3), label = T)) +
  labs(x = "Month", y = "Fecal coliform (log10 MPN/100 ml)")


# Compute summary statistics ----------------------------------------------

# Concentration range in this study compared to recent AA data
summary_this_study <- 
  water %>% 
  filter(target == "ec" & !is.na(location)) %>%
  summarize(q25 = quantile(l10_100ml_cens, probs = 0.25),
            med = median(l10_100ml_cens),
            q75 = quantile(l10_100ml_cens, probs = 0.75)) %>%
  mutate(target = "EC")
summary_AA <- 
  aa %>% 
  filter(station %in% c("Pte Suecia", "Pte La Dehesa", "Pte Bulnes/Panamericana", 
                        "Farfana aguas abajo", "Rinconada de Maipu"),
         date > "2016-01-01") %>%
  summarize(q25 = quantile(l10_fc, probs = .25, na.rm=T),
            med = median(l10_fc, na.rm=T),
            q75 = quantile(l10_fc, probs = 0.75, na.rm=T)) %>%
  mutate(target = "FC")
summary_both <- bind_rows(summary_this_study, summary_AA)


# Write plots/data to file ------------------------------------------------

write_csv(summary_both, here::here("results", "AA_comparison_summary.csv"))
write_tif(plot_hist_data_station, "AA_hist_station.tif")
write_tif(plot_hist_data, "AA_hist_only.tif")
write_tif(plot_hist_and_our_data, "AA_hist_our_data.tif")
write_tif(plot_AA_monthly_trends, "AA_hist_monthly.tif")
write_tif(plot_monthly_trend_2016to18, "AA_hist_monthly_trend_2016to18.tif")
