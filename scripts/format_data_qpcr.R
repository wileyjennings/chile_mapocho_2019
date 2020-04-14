# Doc Info ----------------------------------------------------------------

# Project: Mapocho River, Chile water quality survey
# Author: Wiley Jennings | Boehm Group | Stanford 
# Date: 2 Dec 2019

# Description: Read and format qPCR data. Write .rds file with qPCR data 
# ready to analyze

# Notes: This reads and combines data for all three targets: HF183,
# crAssphage, and norovirus. It includes samples, 1:10 dilutions for testing
# inhibition, and standards, but excludes QC samples and dog and bird
# challenge samples.


# Dependencies ------------------------------------------------------------

required_packages <- c("dplyr", "here")
lapply(required_packages, library, character.only = T)
source(here::here("scripts", "util.R"))


# Read and format ---------------------------------------------------------

hf183_list <- read_s1plus("hf183")
crass_list <- read_s1plus("crass")
noro_list <- read_s1plus("noro")

# Bind all qpcr tibbles together, select data, and format
qpcr <- bind_rows(hf183_list, crass_list, noro_list)  
qpcr <- 
  qpcr %>% 
  select(`Sample Name`, Task, Cт, `Cт SD`, Quantity,  `Quantity SD`, 
         `Ct Threshold`,  HIGHSD, NOAMP, EXPFAIL, date, plate_name, target)
names(qpcr) <- c("name", "task", "ct", "ct_sd", "copy_num", "copy_num_sd", 
                 "ct_threshold", "flag_sd_high", "flag_no_amp", "flag_exp_fail", 
                 "date", "plate_name", "target")
qpcr$ct <- suppressWarnings(as.numeric(qpcr$ct))
qpcr <- qpcr %>% filter(!is.na(task))  # remove rows without data
qpcr <-  # give name to standards
  qpcr %>%
  mutate(name = ifelse(task == "STANDARD", paste0("s", copy_num), 
                       name))

# Add dilution variable, indicating x-fold dilution of an extract: name 
# strings include ":10" if 10-fold dilution. 
# Will produce NA if ":10" not found. Thus, NA indicates dilution factor of 1.
qpcr$diln <- 
  gsub(pattern = ":", replacement = "", 
       str_extract(qpcr$name, ":(\\d)+")) %>%
  as.numeric()
qpcr <- 
  qpcr %>% 
  mutate(diln = ifelse(is.na(diln) & task == "UNKNOWN", 1, diln))

# Remove '1:10' from sample name strings, now that diln variable added
qpcr$name <- gsub(" 1:10", "", qpcr$name)
qpcr$name <- gsub("1:10 ", "", qpcr$name)
qpcr$name <- gsub(" :10", "", qpcr$name)


# Discard standards and samples that will not be used ---------------------

# Discard three samples that were assayed twice for HF183 (s43, s67, s68).
# These samples run on 20191004 should be discarded due to high sd.
# Also discard sewage HF183 since it was run again on the same plate as 
# crAssphage HF183, to increase comparability.
qpcr <- 
  qpcr %>%
  filter(!(target == "hf183" & 
             date == "20191004" & 
             name %in% c("s43", "s67", "s68", "Sewage")))

# Discard samples from CDC UF decay study.
qpcr <- qpcr[-grep(x = qpcr$name, pattern = "Ext ID"), ]

# Discard NTCs and blanks.
qpcr <- qpcr %>% filter(task != "NTC")
qpcr <- qpcr[-grep(x = qpcr$name, pattern = "Blk"), ]
qpcr <- qpcr[-grep(x = qpcr$name, pattern = "Blank"), ]

# Discard dog and bird challenge samples. Note: capitalization was inconsistent.
qpcr <- qpcr[-grep(x = qpcr$name, pattern = "ird"), ]
qpcr <- qpcr[-grep(x = qpcr$name, pattern = "og"), ]

# Discard noro in sewage measurement because this extract was freeze thawed
qpcr <- qpcr %>% filter(!(name == "sewage" & target == "noro"))

# Remove final crass plate run for sewage because standard curve bad and
# sewage estimated previously on another plate.
qpcr <- qpcr %>% filter(!grepl("sewage_20191213", plate_name))

# Split standards from unknowns for writing to separate .rds file. 
# Necessary before estimating summary ct statistics of unknowns.
standards <- qpcr %>% filter(task == "STANDARD")
samples <- qpcr %>% filter(task != "STANDARD")

# Compute mean replicate ct values and compute SDs of unknowns. 
# Discard SDs computed by Step One Plus software.
# Discard replicates now that computed summary values using `distinct()`.
samples <- 
  samples %>%
  group_by(plate_name, target, name, diln) %>%
  summarize(ct_mean = mean(ct, na.rm = T),
            ct_sd = sd(ct, na.rm = T)) %>%
  left_join(x = samples %>% 
              select(-ct_sd) %>% 
              distinct(plate_name, target, name, diln, .keep_all = T), 
            y = ., 
            by = c("plate_name", "target", "name", "diln")) 

# Write data
saveRDS(standards, here::here("data", "processed", "standards.rds"))
saveRDS(samples, here::here("data", "processed", "samples.rds"))
