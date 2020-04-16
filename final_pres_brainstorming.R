
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


# Model 
# compare GAMM and GLMM?


# Results and conclusions
# Did not find strong spatial effect, but did find noticeable temporal effects
