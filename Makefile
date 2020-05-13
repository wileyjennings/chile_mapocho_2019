.PHONY: all
all: format stat_plot

##########################################
#                                        #
#               FORMAT                   #
#                                        #
##########################################

# Set variables: processed data to write
lab = data/processed/lab.rds
prcp = data/processed/prcp.rds
qpcr = data/processed/qpcr.rds
samp = data/processed/samples.rds results/inhibition_summary.csv
water = data/processed/water.rds

# Set variables: scripts
util = scripts/util.R

# Formats data req'd for analysis and writes the data to .rds files
format: $(lab) $(prcp) $(qpcr) $(water)

# Produces prereqs of 'format' rule
$(lab): data/raw/Mapocho_sampling_data.xlsx
	Rscript scripts/01_format_data_lab.R
$(prcp): data/raw/Precip_Santiago_330019_2019_Agua6Horas.csv
	Rscript scripts/01_format_data_prcp.R
# For simplicity, this omits the many qPCR results file dependencies, since\
these files should never be altered
$(qpcr): $(util) 
	Rscript scripts/01_format_data_qpcr.R
$(water): $(lab) $(prcp) $(samp) $(util)
	Rscript scripts/04_format_data_water.R


##########################################
#                                        #
#         QPCR STNDS & CONCS             #
#                                        #
##########################################

# Need to handle issue of multiple targets here
stnd = data/processed/standards_mixed_summary.rds \
data/processed/standards_mixed.rds
out_stnd_eqns = $(stnd) results/standards_summary.csv figures/standards.tif

# Write stnd eqns .rds files as well as summary table and figure
stnd_eqns: $(out_stnd_eqns)

# Write sample concentration .rds and inhibition summary .csv
estim_conc: $(samp)

$(out_stnd_eqns): $(qpcr) $(util)
	Rscript scripts/02_estim_qpcr_stnd_eqns.R
$(samp): $(qpcr) $(stnd)
	Rscript scripts/03_estim_qpcr_rxn_conc.R


##########################################
#                                        #
#        STATS, TABLES, & FIGURES        #
#                                        #
##########################################

stat_plot_scripts = $(wildcard scripts/05*.R)
.PHONY: stat_plot
stat_plot: $(water) $(util)
	$(foreach var,$(stat_plot_scripts),Rscript $(var);)


##########################################
#                                        #
#                  CLEAN                 #
#                                        #
##########################################

.PHONY: clean, clean_data, clean_figures, clean_results

clean: clean_data clean_figures clean_results

clean_data:
	rm -f data/processed/*.rds

clean_figures:
	rm -f figures/*.tif figures/*.png

clean_results:
	rm -f results/*.csv

