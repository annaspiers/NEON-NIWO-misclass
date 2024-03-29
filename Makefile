.PHONY : help clean figures
.DEFAULT_GOAL := figures



# Figures ===========================================================================

## figures : makes figures from main text
figures: figures/thetaconfusion.png \	
figures/site-level-ranefs.png \	
figures/spec-level-ranefs.png \	
figures/occthroughtime.png \		
figures/CIwidthcomparison.png \	
figures/comparedensities.png \	
figures/thetadifference.png \	



# Cleaning Data ====================================================================

## Download and clean NEON data
data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds: carabid_data_clean.R
	R CMD BATCH carabid_data_clean.R



# Analyses ==========================================================================

## Simulation-based callibration for priors

## Full joint occupancy-classification model output
output/full_jm.rds output/full_jmsumm.rds output/full_theta_summ.rds: \
full_dyn_occ_misclass.R source/jags_input.R full_dyn_occ_misclass_JAGS.txt \
data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds
	R CMD BATCH full_dyn_occ_misclass.R
	
	
## Figures for full joint occupancy-classification model 
figures/thetaconfusion.png: create_figs.R source/jags_input.R \
data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds
	R CMD BATCH create_figs.R

figures/site-level-ranefs.png figures/spec-level-ranefs.png: \
create_figs.R output/full_jm.rds 
	R CMD BATCH create_figs.R
	
figures/occthroughtime.png: \
create_figs.R output/full_jm.rds 
	R CMD BATCH create_figs.R
	

## Reduced classification model output
output/reduced_jm.rds output/reduced_jmsumm.rds output/reduced_theta_summ.rds: \
reduced_misclass.R reduced_misclass_JAGS.txt source/jags_input.R \
data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds 
	R CMD BATCH reduced_misclass.R
	
	
## Figures comparing full and reduced models
figures/CIwidthcomparison.png: create_figs.R output/reduced_theta_summ.rds output/full_theta_summ.rds \
source/jags_input.R data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds
	R CMD BATCH create_figs.R
	
figures/comparedensities.png: create_figs.R output/full_jm.rds output/reduced_jm.rds
	R CMD BATCH create_figs.R

figures/thetadifference.png: create_figs.R source/jags_input.R
 data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds
	R CMD BATCH create_figs.R
	

## Validation results 
output/val_full_jm.rds: \
validation.R source/jags_input.R full_dyn_occ_misclass_JAGS.txt \
data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds
	R CMD BATCH validation.R
	
	
## Simulation output
output/full_jm.rds output/full_jmsumm.rds output/full_theta_summ.rds: \
full_dyn_occ_misclass.R source/jags_input.R full_dyn_occ_misclass_JAGS.txt \
data/all_paratax_df.rds data/pinned_df.rds data/expert_df.rds
	R CMD BATCH full_dyn_occ_misclass.R



# Makefile Cleaning ====================================================================================
## clean : removes all files generated by this makefile. Individual steps can be cleaned by using remove-figures 
## which deletes all figures and remove-R-outputs which removes saved R-objects and *.Rout files.

clean: remove-figures remove-R-outputs

remove-figures:
	rm -f figures/*.png

remove-R-outputs:
	rm -f *.Rout
	rm -f Rplots.pdf
	rm -f output/*.rds
	rm -f data/*_df.rds

