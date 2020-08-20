
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
  start_ <- ifelse(.target == "hf183", 113, 118)
  plate_names <- mapply(FUN = substr, x = paths, start = start_, 
                        stop = nchar(paths)-4)
  tib_list <- map2(.x = tib_list, .y = plate_names, ~ mutate(.x, plate_name = .y)) 
  tib_list <- map(.x = tib_list, 
                  ~ mutate(.x, target = ifelse(
                    `Target Name` == "crAss 064", "crass", ifelse(
                      `Target Name` == "HF183", "hf183", "noro")))) 
  
  tib_list
}

# Helper to tidy up model results.
# Given:
# - .data: a tibble containing .model_col, a (list) col of (linear) models.
# Returns:
# - .data with three new list cols with model summaries and fitted values.

tidy_model <- function(.data, .model_col) {
  .model_col <- enquo(.model_col)
  # Mixed model throws annoying class coercion warning. Could be resolved
  # by using broom.mixed functions.
  .data %>%
    mutate(tidied = suppressWarnings(map(!!.model_col, tidy)),
           glanced = map(!!.model_col, glance),
           augmented = map(!!.model_col, augment))
}

# Estimate qPCR standard curves.
# Given:
# - .data: a data frame of qPCR standard measurements containing cols:
# -- `target`, `ct`, `copy_num`, and `plate_name`.
# Returns:
# - A tibble containing a list column of random effects models, one model
# -- for each qPCR target; a list col of random effects; and list cols of 
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
    mutate(
      mixed = map(
        data, ~lmer(ct ~ log10(copy_num) + (1|plate_name), data = .x)),
      re = map(mixed, ranef),
      re = map(re, function(x) as_tibble(as.list(x))))

  # Return nicely formatted mixed effects model results
  tidy_model(model, mixed)
}

# Summarize standard equation parameters.
# Given:
# - A tibble containing a column of tidied models (e.g., from `tidy_model()`)
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
  # Return summary with qPCR efficiency col and conditional R2
  std_summary %>%
    mutate(
      effic = 10^(-1/slope_l10cn) - 1,
      r2_cond = map(standards_mixed$mixed, performance::r2),
      r2_cond = sapply(r2_cond, "[[", 1))
}

# Estimate lower limit of detection using logistic regression to estimate
# the standard concentration that amplified with probability `.prob`
# Given:
# - .data: data frame of standards
# - .target: the qPCR target, a string
# - .prob: numeric [0, 1], the lod corresponds to the lowest copy num with at 
# -- least .prob probability of amplifying.
# Returns:
# - Numeric LOD, in units of log10(copy number)
lod_logit <- function(.data, .target, .prob) {
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
  
  # Return lod in units of copy number
  grid %>%
    filter(abs(.prob - pred) == min(abs(.prob - pred))) %>%
    pull(l10_cn)
}

# Estimate coefficient of variation for log-transformed qPCR data,
# derived by Forootan 2017
# Given:
# - ct_sd: a vector of Ct standard deviations
# - eff: the qpcr efficiency
# Returns:
# - a vector of coefficient of variations with same length as ct_sd
cv_qpcr <- function(ct_sd, eff = 1) {
  power <- ct_sd^2 * log(1+eff)
  sqrt((1 + eff)^power - 1)
}

# Estimate lower limit of quantification, namely the lowest concentration that 
# can be precisely quantified, based on the coefficient of variation.
# Pooled standard deviation across plates is used in the CV calculation;
# CV calculation is based on Cq's (see function: cv_qpcr)
# Given:
# - .data: data frame of standards
# - .target: the qPCR target, a string
# - .eff: the qpcr efficiency for that target, a number
# Returns:
# - a list containing 3 elms: 
# -- 1. data frame containing the stnd copy_num, name, and CV
# -- 2. Predicted CV as a function of log10 copy num, interpolated using spline
# -- 3. A plot object to plot CV as function of copy_num, with smooth line
lloq_pool <- function(.data, .target, .eff) {
  data <- 
    .data %>% 
    filter(target == .target, !is.na(ct))  # removes wells that did not amplify
  
  # Count the number of standards run at each conc on each plate for pooled sd
  # estimation
  count_std <- 
    data %>%
    count(plate_name, copy_num)
  
  # Compute sd for each plate
  sd_plate <- 
    data %>%
    group_by(copy_num, plate_name) %>%
    summarize(sd = sd(ct)) %>%
    left_join(., count_std, by = c("plate_name", "copy_num"))
  
  # Compute pooled sd across plates and use that for CV calculation
  cv <- 
    sd_plate %>%
    filter(!is.na(sd)) %>%
    mutate(sd_n = sd*n) %>%
    group_by(copy_num) %>%
    summarize(sd_n_sum = sum(sd_n),
              n_sum = sum(n),
              sd_pool = sd_n_sum / n_sum) %>%
    mutate(cv = cv_qpcr(ct_sd = sd_pool, eff = .eff)) %>%
    select(copy_num, cv) %>%
    mutate(target = .target)
  
  # Interpolate CV using spline fit over grid of log10 copy_num values
  spline_fit <- smooth.spline(
    x = log10(cv$copy_num), y = cv$cv, spar = 0.3)
  spline_pred <- 
    predict(spline_fit, 
            x = seq(min(log10(cv$copy_num)),
                    max(log10(cv$copy_num)), 
                    0.1)) %>%
    as_tibble()
  names(spline_pred) <- c("copy_num_l10", "cv_spline_fit")
  
  # Compute derivative of spline fit for LOQ determination
  spline_pred$deriv <- 
    predict(spline_fit,
            x = seq(min(log10(cv$copy_num)),
                    max(log10(cv$copy_num)), 
                    0.1),
            deriv = 1)$y
  
  # Plot CV against copy_num, with smoothing spline
  loq_plot <- 
    cv %>%
    ggplot(., aes(x = log10(copy_num), y = cv)) +
    geom_point() +
    geom_line(data = spline_pred, aes(x = copy_num_l10, y = cv_spline_fit)) +
    geom_hline(aes(yintercept = 0.35), lty = 2, color = "red") +
    labs(x = "Copy number per rxn, log10", 
         y = "Coefficient of variation",
         title = .target) +
    theme_bw()
  
  # LOQ is the lowest concentration with a CV < 0.35 and a negative slope
  loq_l10 <- 
    spline_pred %>%
    filter(cv_spline_fit < 0.35, deriv < 0) %>%
    summarize(loq_l10 = min(copy_num_l10)) %>%
    pull(loq_l10)
  
  # Return a list of these outputs
  list(
    cv = cv,
    spline_pred = spline_pred,
    loq_plot = loq_plot,
    loq_l10 = loq_l10)
}


# Calculate microbial concentration per 100 ml water from qPCR rxn 
# concentrations
# Given:
# - .data: data frame (`water`) containing qpcr rxn concentrations and lab 
# -- processing values ...
# - .l10_conc_100ml: column name of .data containing l10 concentration per
# -- 100 ml; must already exist in .data.
# - .l10_conc: l10 concentration in qpcr reaction
# Returns:
# - The data frame .data with the column .l10_conc_100ml filled in for qpcr
# -- method and also with .l10_conc column removed, since we are no longer 
# -- interested in the qpcr reaction concentrations.
calc_qpcr_water <- function(.data, .l10_conc_100ml, .l10_conc) {
  # Tidy evaluation
  .l10_conc_100ml <- enquo(.l10_conc_100ml)
  .l10_conc_100ml_name <- quo_name(.l10_conc_100ml)
  .l10_conc <- enquo(.l10_conc)
  
  # (cn/ml water) = 
  #   (cn/rxn)*(1 rxn/5 ul extract)*(100 ul extract/1 extraction)*
  #   (1 extraction/X ml water filtered)
  # (cn/100 ml water) = (cn/ml water)*100
  
  vol_extract <- 100  # 100 ul nucleic acid extract
  vol_extract_qpcr <- 5  # 5 ul nucleic acid used as template in qpcr rxns
  # .l10_conc is concentration in cn/rxn = cn/5 ul extract
  
  .data %>%
    mutate(
      !!.l10_conc_100ml_name := ifelse(
        method == "qpcr",
        log10(
          (10^!!.l10_conc)*(vol_extract/vol_extract_qpcr)*(100/vol_mce)),
        !!.l10_conc_100ml)) %>%
    select(-!!.l10_conc)
}

# Format model residuals as time series.
# Given .data containing a list column of models, a .target column, and a list
# col of .resids residuals
# Returns a vector of time series formatted residuals.
ts_resids <- function(.data, .target, .resids) {
  .resids <- enquo(.resids)
  .data %>% 
    filter(target == .target) %>% 
    pull(!!.resids) %>% 
    unlist(.) %>%
    ts(.)
}


# Plotting & writing functions ------------------------------------------------

# Write a .tif image to file.
# Given:
# - .fig: A plot object
# - .file_name: the name of the file (not path); a string.
# Returns:
# - Nothing.
# Side effect:
# - Writes a high res .tif file to figures/
write_tif <- function(.fig, .file_name, .orient = "horiz") {
  if(!(.orient %in% c("horiz", "vert"))) {
    stop("`orient` must be `horiz` or `vert`")}
  if(.orient == "horiz") {
    .width = 6
    .height = 4
  } else {
    .width = 4
    .height = 6
  }
  tiff(here::here("figures", .file_name), width = .width, height = .height, 
       units = "in", res = 600)
  print(.fig)
  dev.off()
}

# Plot facetted environmental vars
# Given: a ggobj with defined x and y aesthetics, as well as
# geometries (e.g., geom_point)
# Returns: a ggobj prettied up, including colors, text size,
# and desired facet format.
facet_env <- function(.ggobj, .nrow, .num_targ) {
  color_order <- c(1, 4, 2, 3, 5, 6)[1:.num_targ]
  .ggobj +
    facet_wrap(~meas_class, nrow = .nrow, scales = "free_y", 
               strip.position = "left") +
    scale_color_manual(values = viridis::viridis_pal(
      end = 0.90, option = "C")(.num_targ)[color_order]) +
    scale_fill_manual(values = c("gray100", "gray 25")) +
    guides(
      color = guide_legend(title = "", label.hjust = 0.1, 
                           keywidth = unit(0.01, "cm")),
      fill = guide_legend(
        title = "", label.hjust = 0.1, keywidth = unit(0.01, "cm"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = c(1.11, 0.7),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.3, 2, 0, 0), "cm"),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 9),
          legend.text = element_text(size = 8))
}

# Write a table to .csv file, using the name of the object as the file name.
write_result <- function(.obj) {
  .obj_name <- quo_name(enquo(.obj))
  write_csv(.obj, here::here("results", paste0(.obj_name, ".csv")))
}
