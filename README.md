# Analysis of Mapocho River (Santiago, Chile) microbial water quality data.

Data were collected for research purposes in 2019. The PI is Alexandria Boehm, Dept. of Civil and Environmental Engineering, Stanford University. All analysis was performed by Wiley Jennings.

These scripts may be of interest to others performing analysis of environmental sampling data, for which censoring limits (e.g., LOD, LLOQ), spatial clustering, and serial correlation are common statistical issues.

This repository is structured as a full reproducible R project. It contains the raw data, as well as directory structure and R scripts, to reproduce all analyses. 

 Scripts can be run individually, in interactive mode or from command line (scripts are written such that they do not require any input arguments). In this case, they should be run in order, from scripts starting `01` to `05`.

To reproduce the entire analysis or groups of scripts with a single command, execute the `Makefile`. To reproduce all analyses, at the command line simply type `make`. To execute the subset of scripts that formats all data but does not produce statistical analysis or plots, type `make format`. To remove all files produced from running these scripts, type `make clean`. Other clean rules are also written in `Makefile`.

Dependencies are managed with packrat, so that required package versions are included in this repo.
