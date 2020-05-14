# Doc Info ----------------------------------------------------------------

# Title: Pseudo makefile
# Project: Chile Mapocho survey
# Author: Wiley Jennings | Boehm Group | Stanford
# Date: 23 Apr 2020

# Description: pseudo make function to run groups of scripts from command line. 
# You can run:
# --all: all the scripts to reproduce the entire analysis including figures, 
# --format: just the scripts up through level '04' required to produce the `water`
# data frame that serves as the basis for analysis and plotting, or
# --analysis: the level `05 scripts` to produce analysis, assuming all scripts up 
# through `04` have been run.

make <- function() {
  args <- commandArgs(trailingOnly = T)
  args_allowed <- c("--all", "--format", "--analysis")
  if(length(args) != 1 | !(args[1] %in% args_allowed)) {
    stop("\nUsage:  Rscript main.R --<arg>\n  where <arg> is 'all', 'format', or 'analysis'.")
    }
  
  scripts_format <- Sys.glob(here::here("scripts", "0[1-4]*"))
  scripts_analysis <- Sys.glob(here::here("scripts", "05*"))
  
  if(args[1] == "--format") {lapply(scripts_format, source)} 
  else if(args[1] == "--analysis") {lapply(scripts_analysis, source)}
  else{lapply(c(scripts_format, scripts_analysis), source)}
}

make()
         