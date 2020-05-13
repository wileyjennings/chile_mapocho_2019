# Analysis of Mapocho River (Santiago, Chile) microbial water quality data.

Data were collected for research purposes in 2019. The PI is Alexandria Boehm, Dept. of Civil and Environmental Engineering, Stanford University. All analysis was performed by Wiley Jennings.

These scripts may be of interest to others performing analysis of environmental sampling data, for which censoring limits (e.g., LOD, LLOQ), spatial clustering, and serial correlation are common statistical issues.

This repository is structured as a full reproducible R project. It contains the raw data, as well as directory structure and R scripts, to reproduce all analyses. 

On a UNIX system, the entire analysis or groups of scripts can be run with a single command using the `Makefile`. Specifically, to reproduce all analyses, at the command line simply type `make`. To execute the subset of scripts that formats all data but does not produce statistical analysis or plots, type `make format`. To remove all output files produced by these scripts, type `make clean`. Other clean rules are also written in `Makefile`.

On other OS, you will have to run scripts individually. If you run them from the command line, simply use `Rscript <script_name>` (scripts are written such that they do not require any input arguments). Scripts must be run in order, from scripts starting `01` to `05`.

Dependencies are managed with packrat, so that required package versions are included in this repo.
