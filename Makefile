##########################################
#                                        #
#               FORMAT                   #
#                                        #
##########################################

# Set variables: processed data
lab = data/processed/lab.rds
prcp = data/processed/prcp.rds
qpcr = data/processed/qpcr.rds

# Set variables: scripts
util = scripts/util.R

# Formats and writes formatted data to .rds files
format: $(lab) $(prcp) $(qpcr)

$(lab): data/raw/Mapocho_sampling_data.xlsx
	Rscript scripts/01_format_data_lab.R
$(prcp): data/raw/Precip_Santiago_330019_2019_Agua6Horas.csv
	Rscript scripts/01_format_data_prcp.R
# For simplicity, this omits the many qPCR results file dependencies
$(qpcr): $(util) 
	Rscript scripts/01_format_data_qpcr.R


##########################################
#                                        #
#         QPCR STNDS & CONCS             #
#                                        #
##########################################

# Need to handle issue of multiple targets here
out_stnd_eqns = data/processed/standards_mixed_summary.rds \
data/processed/standards_mixed.rds \
results/standards_summary.csv \
figures/standards.tif

stnd_eqns: $(out_stnd_eqns)

$(out_stnd_eqns): $(qpcr) $(util)
	Rscript scripts/02_estim_qpcr_stnd_eqns.R





##########################################
#                                        #
#        STATS, TABLES, & FIGURES        #
#                                        #
##########################################
