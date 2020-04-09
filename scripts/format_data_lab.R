# Doc Info ----------------------------------------------------------------

# Project: Mapocho River, Chile water quality survey
# Author: Wiley Jennings | Boehm Group | Stanford 
# Date: 2 Dec 2019

# Description: Read and format non-qPCR lab data.

# Notes:
# Lab data: sample names and ids, sample meta-data, lab book keeping, and
# - non-qPCR water quality measurements made in lab
# - the lab data file is called 'sampling data'
# - lab$sample_id is equivalent to qpcr$name
# - hourly sampling campaign has FIB but not qPCR measurements


# Dependencies ------------------------------------------------------------

required_packages <- c("here", "lubridate", "readxl", "tidyverse")
lapply(required_packages, library, character.only = T)


# Read and format ---------------------------------------------------------

lab <- 
  read_excel(here::here("data", "raw", "Mapocho_sampling_data.xlsx"),
             sheet = "all_samples",
             col_types = c(rep("text", 5),
                           rep("numeric", 12),
                           "text",
                           rep("numeric", 3),
                           "text",
                           rep("numeric", 3),
                           "text",
                           rep("numeric", 2),
                           rep("text", 7)), 
             na = c("NA", ""))

# Remove IDEXX well counts
lab <- lab[, -c( grep("large", names(lab)), grep("small", names(lab)) )] 

names(lab) <- c("campaign", "sample_name", "sample_id",
                "date", "time",  "turb", "cond", "temp",
                "vol_mce", "vol_pc", "tc_conc", "tc_cens",
                "tc_hi", "tc_lo", "ec_conc", "ec_cens", "ec_hi",
                "ec_lo", "ent_conc", "ent_cens", "ent_hi", "ent_lo",
                "notes", "hf183_assayed", "hf183_inhib_test",
                "crass_assayed", "crass_inhib_test", 
                "noro_assayed", "noro_inhib_test")

# Round lab time to nearest half hour. This is so that temporal campaign sample
# times match at 30 minute resolution
lab$datet <- as.POSIXct(paste(lab$date, lab$time), tz = "America/Santiago", 
                        format = "%m/%d/%Y %H:%M")
lab$datet <- round_date(lab$datet, "30 minutes")

# Reshape lab to long format for processing microbial data with dplyr
lab <-
  lab %>%
  pivot_longer(cols = c(tc_conc:ent_lo), 
               names_to = c("target", ".value"), 
               names_sep = "_")

# Provide more descriptive names for microbial concentration vars
lab <- lab %>% rename(conc_100ml_cens = conc, hi95 = hi, lo95 = lo)

# Log transform data and define censoring indicator variables
lab <- 
  lab %>%
  mutate(l10_turb = log10(turb)) %>%
  select(-turb)
lab <- 
  lab %>%
  mutate(l10_turb = log10(turb),
         l10_100ml_cens = log10(conc_100ml_cens),
         l10_100ml_hi95 = log10(hi95),
         l10_100ml_lo95 = log10(lo95),
         lod_l10 = 1,
         loq_l10 = NA,
         cens = ifelse(cens == "left", "blod", 
                       ifelse(cens == "no", "ncen", "aloq"))) %>%
  select(-turb)

# Remove observations where no samples was taken
lab <- lab %>% filter(!grepl("NO SAMPLE", sample_name))