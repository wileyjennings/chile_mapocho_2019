
# Requirements ------------------------------------------------------------

required_packages <- c("here", "tidyverse")
lapply(required_packages, library, character.only = T)


# Reads qPCR data files in default format produced by Step One Plus qPCR platform.

# Given:
# - target: the qPCR target. Can be "HF183", "crass", or "noro"
# --  a string
# Returns:
# - A list of tibbles, where each element is the results from one instrument
# - run. The plate name and date have been appended to each tibble element to
# - distinguish them.

# Requirements
# - qpcr .xls files which are titled as:
# -- Begins with target name.
# -- Contains date formatted as YYYYMMDD.

read_s1plus <- function(target_) {
  file_name <- paste0(target_, "*.xls")
  path <- here::here("data", "raw", "qPCR", target_, "results", 
                     paste0(target_,"*.xls"))
  paths <- Sys.glob(path)
  
  # Read results sheets from qPCR workbook files
  tib_list <- map(paths, read_excel, sheet = "Results", skip = 7)
  
  # Add date to each plate run
  dates <- str_extract(paths, "2019(\\d)+")
  tib_list <- map2(.x = tib_list, .y = dates, ~ mutate(.x, date = .y)) 
  
  # Add plate name and target to each plate run 
  if (target == "HF183") {
    start_ <- 99
  } else if (target %in% c("crass", "noro")) {
    start_ <- 104
  } else {stop("Target must be 'HF183', 'crass', or 'noro'")}
  plate_names <- mapply(FUN = substr, x = paths, start = start_, 
                        stop = nchar(paths)-4)
  tib_list <- map2(.x = tib_list, .y = plate_names, ~ mutate(.x, plate_name = .y)) 
  tib_list <- map(.x = tib_list, ~ mutate(.x, target = target_)) 
  
  tib_list
}


