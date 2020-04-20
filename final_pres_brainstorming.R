
# Brainstorming for forecasting final presentation

# Premise


# What were our questions
# Hypothesized strong spatial/microhabitat effects
# Predict when is peak abundance for each species in 2019. Assess based on timing and magnitude of peak. Area under curve for total seasonal abundance
# Here, we focus on two species. Carabus taedatus (more in tundra) and Cymindis unicolor (more in fores) are the only two species present in both nlcdClass


# EDA
## Beetles

#   Identification - we already showed this at one of the earlier classes, so we might not need to show this now
# 
#   Species selection - 2/7 species that were most abundant
#       ?info/natural history about the species we chose
# 
#   Abundance over time, per trap - awesome figure Anna made
# 
#   ?Map of beetle traps
# 
# 
## Predictors
# 
#   fixed effects:
#     day of year
#   nlcd class, elevation - included with beetle data product
# 
#   canopy height - related to nlcd class and elevation (treeline)
#       raster plot
#       buffer radius choice
#       show scatterplots of correlations?
#     
#   leaf area index
#       raster plot
#       buffer radius choice
# 
#   ?degree days, precipitation - short term and long term
#   ?snowmelt date
# 
# 
#   desired but unused: 
#   soil moisture - not colocated with traps, ?map
#   woody debris - not colocated with traps, ?map
#   vegetation - difficult/unclear to summarize at trap level

# Beetles...
# 
# The literature supported certain explanatory variables, but we didn't some of these to explain much variance in our models.
# 
#  Predictor variables...
# How did we calculate LAI (and other AOP data products) for each trap?
    # How did we decide on what size buffer to use for these calculations?

# ANNA'S STORYTELLING BRAINSTORM
# Out of these accurately id abundant species, many were present in only one
# habitat, only these two were present in both, wanted to explore more. through
# the sampling - neon found that one of the plots had low abundance (4) so
# switched out for plot 13, can we predict the abundance for plot 13 in 2018 for
# these two species??? given that we were intereste in these two species, we dug
# into their biology....natural history of both species chosen....for these
# reasons, we chose these explanatory variables to include....this is what the
# process was like to summarize each predictor variable for our model (e.g.
# defining buffer for canopy height)...heres a visual of the two species through
# time...this is what our model was...questions: what predictors are most
# important for explaining abundance of these species...what was our prediction
# of the abundance, how accurate, how well does it match the raw data?...


# Model 
# compare GAMM and GLMM?


# Results and conclusions
# Did not find strong spatial effect, but did find noticeable temporal effects
