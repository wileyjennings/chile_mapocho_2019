
# Requirements ------------------------------------------------------------

required_packages <- c("broom", "dplyr", "glue", "here", "purrr", "readxl", 
                       "stringr", "tidyr")
lapply(required_packages, library, character.only = T)


# Functions ---------------------------------------------------------------

# Reads qPCR data files in default format produced by Step One Plus qPCR platform.
# Given:
# - target: the qPCR target. Can be "HF183", "crAssphage", or "noroGII"
# --  a string
# Returns:
# - A list of tibbles, where each element is the results from one instrument
# - run. The plate name and date have been appended to each tibble element to
# - distinguish them.
# Requirements
# - qpcr .xls files which are titled as:
# -- Begins with target name.
# -- Contains date formatted as YYYYMMDD.

read_s1plus <- function(.target) {
  if (.target == "hf183") {
    path_target <- "HF183"
  } else if (.target == "crass") {
    path_target <- "crAssphage"
  } else if (.target == "noro") {
    path_target <- "noroGII"
  } else {stop("Target must be 'hf183', 'crass', or 'noro'")}
  file_name <- paste0(path_target, "*.xls")
  path <- here::here("data", "raw", "qPCR", path_target, "results", 
                     paste0(path_target,"*.xls"))
  paths <- Sys.glob(path)
  
  # Read results sheets from qPCR workbook files
  tib_list <- map(paths, read_excel, sheet = "Results", skip = 7)
  
  # Check that `Target Name` is expected
  get_targets <- function(df) {
    df %>%
      distinct(`Target Name`) %>%
      filter(!is.na(`Target Name`)) %>%
      pull()
  }
  target_names <- map(tib_list, get_targets) %>% unlist()
  if(!all(target_names %in% c("crAss 064", "HF183", "HNV GII_ORF1/ORF2_SB"))) {
    stop("`Unexpected `Target Name`!")
  }
  
  # Add date to each plate run
  dates <- str_extract(paths, "2019(\\d)+")
  tib_list <- map2(.x = tib_list, .y = dates, ~ mutate(.x, date = .y)) 
  
  # Add plate name and target to each plate run 
  start_ <- ifelse(.target == "hf183", 99, 104)
  plate_names <- mapply(FUN = substr, x = paths, start = start_, 
                        stop = nchar(paths)-4)
  tib_list <- map2(.x = tib_list, .y = plate_names, ~ mutate(.x, plate_name = .y)) 
  tib_list <- map(.x = tib_list, 
                  ~ mutate(.x, target = ifelse(
                    `Target Name` == "crAss 064", "crass", ifelse(
                      `Target Name` == "HF183", "hf183", "noro")))) 
  
  tib_list
}

# Tidy up model results.
# Given:
# - .data: a tibble containing .model_col, a (list) col of (linear) models.
# Returns:
# - .data with three new list cols with model summaries and fitted values.

tidy_model <- function(.data, .model_col) {
  .model_col <- enquo(.model_col)
  .data %>%
    mutate(tidied = map(!!.model_col, tidy),
           glanced = map(!!.model_col, glance),
           augmented = map(!!.model_col, augment))
}

# Estimate qPCR standard curves.
# Given:
# - .data: a data frame of qPCR standard measurements containing cols:
# -- `target`, `ct`, `copy_num`, and `plate_name`.
# Returns:
# - A tibble containing a list column of random effects models, one model
# -- for each qPCR target, a list col of random effects, and list cols of 
# -- nicely tidied model results.

estim_std_mix <- function(.data) {
  if(!all(c("target", "ct", "copy_num", "plate_name") %in% names(.data))) {
    stop("`.data` must contain `target`, `ct`, `copy_num`, and `plate_name` 
         columns")
  }
  
  # Warn if few clusters for estimating random effects
  few_clusters <- 
    .data %>% 
    group_by(target) %>% 
    distinct(plate_name) %>% 
    summarize(n = n()) %>% 
    filter(n < 6) 
  cluster_num <- few_clusters %>% pull(n)
  cluster_target <- few_clusters %>% pull(target)
  if (length(cluster_target) > 0) {
    warning(glue(
      "{cluster_target} was run on only {cluster_num} plates. \nAre you sure you want to use a mixed model approach?\n"
    ))
  } 
  
  # Estimate random effects model
  model <- 
    .data %>%
    nest(data = -target) %>%
    mutate(mixed = map(data, ~lmer(ct ~ log10(copy_num) + (1|plate_name), 
                                   data = .x)),
           re = map(mixed, ranef),
           re = map(re, function(x) as_tibble(as.list(x))))
  
  # Return nicely formatted mixed effects model results
  tidy_model(model, mixed)
}

# Summarize standard equation parameters.
# Given:
# - A tibble containing a column of tidied models (e.g., from `tidy_mode()`)
# Returns:
# - A tibble containing important standard equation parameters.
summarize_std <- function(.data) {
  if(!("tidied" %in% names(.data))) {
    stop("`.data` must contain `tidied` column")
  }
  
  std_summary <- 
    .data %>%
    unnest(tidied) %>%
    pivot_wider(., id_cols = target, names_from = term, 
                values_from = c(estimate))
  names(std_summary) <- c("target", "intercept", "slope_l10cn","intercept_sd", 
                          "resid_sd")
  # Return summary with qPCR efficiency col
  std_summary %>%
    mutate(effic = 10^(-1/slope_l10cn) - 1)
}

# Estimate lower limit of quantification using logistic regression to estimate
# the standard concentration that amplified with probability `.prob`
# Given:
# - .data: data frame of standards
# - .target: the qPCR target, a string
# - .prob: numeric [0, 1], the lloq corresponds to the lowest copy num with at 
# -- least .prob probability of amplifying.
# Returns:
# - Numeric LLOQ, in units of copy number
lloq_logit <- function(.data, .target, .prob) {
  # Add presence absence col
  .data <- 
    .data %>%
    filter(target == .target) %>%
    mutate(amp_bin = ifelse(is.na(ct), 0, 1),
           l10_cn = log10(copy_num)) 
  
  logist_mod <- glm(amp_bin ~ l10_cn, data = .data, family = "binomial")
  
  # Predict detection probability on grid of copy numbers
  grid <- tibble(cn = c(1:1000), l10_cn = log10(cn))
  grid$pred <- predict(logist_mod, newdata = grid, type = "response")
  
  # Return lloq in units of copy number
  grid %>%
    filter(abs(.prob - pred) == min(abs(.prob - pred))) %>%
    pull(cn)
}

