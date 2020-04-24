# Analysis of Mapocho River (Santiago, Chile) microbial water quality data.

Data were collected for research purposes in 2019. The PI is Alexandria Boehm, Dept. of Civil and Environmental Engineering, Stanford University. All analysis was performed by Wiley Jennings.

These scripts may be of interest to others performing analysis of environmental sampling data, for which censoring limits (e.g., LOD, LLOQ), spatial clustering, and serial correlation are common statistical issues.

This repository is structured as a full reproducible R project. It contains the raw data, as well as directory structure and R scripts, to reproduce all analyses. 

 Scripts can be run individually, in interactive mode or from command line (scripts are written such that they do not require any input arguments). In this case, they should be run in order, from scripts starting `01` to `05`.

To reproduce the entire analysis immediately, scripts can also be run using the main file, which allows all formatting scripts to be run with one command, and then analysis scripts to be run with another, or for all scripts to be run in order. Simply fork the repository, and then from command line in the project directory, type
  Rscript scripts/main.R
which will prompt an error providing proper usage.

Dependencies are managed with packrat, so that packages do not need to be installed, and the proper package version is included in the repo. 
