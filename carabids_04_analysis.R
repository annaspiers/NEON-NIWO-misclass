# We try different models of predicting carabid abundance and composition

library(lme4)      #max lik multilevel: lmer(), glmer() etc
library(arm)       #for se.ranef()
library(ggplot2)
library(gridExtra) #arranging multiple plots
library(dplyr)
library(mgcv)       #for gams
library(brms)
library(mctest)     #omcdiag()
library(gratia)
library(jagsUI)

all7sp_dat <- read.csv("data_derived/model_df_by_species_in_sample.csv") %>%
    mutate(sc_DOY = scale(DOY, center = TRUE, scale = TRUE),
           col_year_fac = as.factor(col_year))

# Check for collinearity of predictors ----------------------------------------------------------

# resource: https://cran.r-project.org/web/packages/olsrr/vignettes/regression_diagnostics.html
# resource: https://datascienceplus.com/multicollinearity-in-r/
# VIF=1 means no variance of predictor. VIFs >4 warrants investigation. VIFs >10 indicate serious multicollinearity
vars_mat <- all7sp_dat %>%
    mutate(nlcdclass_num = as.numeric(nlcdClass)) %>%
    dplyr::select(LAI_1718avg, trap_CHM, DOY, col_year, nlcdclass_num)#, elev )
omcdiag(vars_mat, all7sp_dat$sp_abund)
imcdiag(vars_mat, all7sp_dat$sp_abund, method="VIF") #failed to dect multicollinearity when elevation was remoted
# CHM and nlcdClass are highly correlated, but CHM contains more info, so leave out nlcdClass. nlcdClass also had a large VIF


# Compare GAMM and GLMM ---------------------------------------------------

# For Carabus taedatus, we see higher abundance in the tundra 

# all variables
glmm1_cartae <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY + nlcdClass +
                      (1|plotID:trapID) + (1|col_year:collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus") )
# with squared sc_DOY
glmm2_cartae <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY^2 + nlcdClass +
                      (1|plotID:trapID) + (1|col_year:collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus") )
# with sc_DOY, no col_year grouping
glmm3_cartae <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY + nlcdClass +
                      (1|plotID:trapID) + (1|collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus") )
# with squared sc_DOY, no col_year grouping
glmm4_cartae <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY^2 + nlcdClass +
                      (1|plotID:trapID) + (1|collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus") )
AIC(glmm1_cartae)
AIC(glmm2_cartae)
AIC(glmm3_cartae)
AIC(glmm4_cartae)


# all variables with splines (error when nlcdClass has spline)
gamm1_cartae <- gam(sp_abund ~ s(LAI_1718avg) + s(trap_CHM) + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus"))
# all variables remove splines on predictors except DOY
gamm2_cartae <- gam(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus"))
# remove nlcdClass
gamm3_cartae <- gam(sp_abund ~ LAI_1718avg + trap_CHM  + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus"))
# remove col_year grouping variable w nlcdClass
gamm4_cartae <- gam(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus"))
# remove col_year grouping variable w/o nlcdClass
gamm5_cartae <- gam(sp_abund ~ LAI_1718avg + trap_CHM  + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus"))
AIC(gamm1_cartae) #tied for first
AIC(gamm2_cartae) #tied for first
AIC(gamm3_cartae)
AIC(gamm4_cartae) #tied for first
AIC(gamm5_cartae)

# Now for the second species

# all variables
glmm1_cymuni <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY + nlcdClass +
                      (1|plotID:trapID) + (1|col_year:collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor") )
# with squared sc_DOY
glmm2_cymuni <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY^2 + nlcdClass +
                      (1|plotID:trapID) + (1|col_year:collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor") )
# with sc_DOY, no col_year grouping
glmm3_cymuni <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY + nlcdClass +
                      (1|plotID:trapID) + (1|collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor") )
# with squared sc_DOY, no col_year grouping
glmm4_cymuni <- glmer(sp_abund ~ LAI_1718avg + trap_CHM + sc_DOY + nlcdClass +
                      (1|plotID:trapID) + (1|collectDate),
                  family=poisson(), data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor") )
AIC(glmm1_cymuni)
AIC(glmm2_cymuni)
AIC(glmm3_cymuni)
AIC(glmm4_cymuni) #all equal


# all variables with splines (error when nlcdClass has spline)
gamm1_cymuni <- gam(sp_abund ~ s(LAI_1718avg) + s(trap_CHM) + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor"))
# all variables remove splines on predictors except DOY
gamm2_cymuni <- gam(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor"))
# remove nlcdClass
gamm3_cymuni <- gam(sp_abund ~ LAI_1718avg + trap_CHM  + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor"))
# remove col_year grouping variable w nlcdClass
gamm4_cymuni <- gam(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor"))
# remove col_year grouping variable w/o nlcdClass
gamm5_cymuni <- gam(sp_abund ~ LAI_1718avg + trap_CHM  + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor"))
AIC(gamm1_cymuni) 
AIC(gamm2_cymuni) #tied for first
AIC(gamm3_cymuni)
AIC(gamm4_cymuni) #tied for first
AIC(gamm5_cymuni)

# For both species, gamms performed better than glmms, and models with nlcdClass and only splining DOY performed best


# Compare GAM and GAMM functions
gam_cartae <- gam(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year_fact, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
gamM_cartae <- gamm(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4),
                   random=list(plotID=~1, plot_trap=~1, col_year_fact=~1, collectDate=~1),
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
gam.check(gamm_cartae)
summary(gamm_cartae)
plot(gamm_cartae, pages=1)


gam_cymuni <- gam(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year_fact, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
gamM_cymuni <- gamm(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4),
                   random=list(plotID=~1, plot_trap=~1, col_year_fact=~1, collectDate=~1),
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
gam.check(gamm_cymuni)
summary(gamm_cymuni)
plot(gamm_cymuni, pages=1)


# Add in precip and slope/aspect to best gam ------------------------------

gam_cartae_wo <- gam(sp_abund ~ LAI_1718avg + trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
gam_cartae <- gam(sp_abund ~ s(precip_2weeks) + trap17aspect +trap17slope + LAI_1718avg +
                      trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
gam_cartae_int <- gam(sp_abund ~ s(precip_2weeks) + te(trap17aspect,trap17slope) + LAI_1718avg +
                      trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
# Model with slope/aspect interaction has the lowest REML score and most deviance explained, but BARELY

# Try plot slope/aspect now (trap before)
gam_cartae_int_plot <- gam(sp_abund ~ s(precip_2weeks) + te(plot17aspect,plot17slope) + LAI_1718avg +
                      trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
gam_cartae_precip_doy <- gam(sp_abund ~ te(DOY, precip_2weeks) + te(plot17aspect,plot17slope) + LAI_1718avg +
                      trap_CHM + nlcdClass +
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
stepAIC(gam_cartae)

gam_cymuni_wo <- gam(sp_abund ~  LAI_1718avg +
                      trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
gam_cymuni_add <- gam(sp_abund ~  (precip_2weeks) + trap17aspect + trap17slope + LAI_1718avg +
                      trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
gam_cymuni_int_trap <- gam(sp_abund ~  s(precip_2weeks) + te(trap17aspect,trap17slope) + LAI_1718avg +
                      trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
gam_cymuni_int_plot <- gam(sp_abund ~  s(precip_2weeks) + te(plot17aspect,plot17slope) + LAI_1718avg +
                      trap_CHM + nlcdClass + s(DOY,  bs = "cc", k=4) + 
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
gam_cymuni_plot_doy <- gam(sp_abund ~ te(DOY, precip_2weeks) + te(plot17aspect,plot17slope) +
                      LAI_1718avg + trap_CHM + nlcdClass +
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=poisson(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))


# Take best model for each species ----------------------------------------

# It is important for a given study to determine if there is evidence for a sig- nificant global smooth effect, we recom- mend fitting model GS or GI, including the argument select = TRUE in the gam function. This has the effect of adding an extra penalty to each smooth term, that penalizes functions in the null space of the penalty matrices for each smooth. By doing this, it is possible for mgcv to penalize all model terms to a zero effect, in effect doing variable selection (Marra & Wood, 2011). When select=TRUE, the significance of the global term can be found by looking at the significance of the term in summary.gam(model). Note that this can significantly increase the amount of time it takes to fit a model for data sets with a large number of penalty terms (such as model GI when the number of groups is high


# Carabus taedatus
# What grouping variables explain the most variance? 
gam_cartae_re <- gam(sp_abund ~ s(collectDate, bs="re") + 
                      s(col_year_fac, bs="re") + 
                      s(plot_trap, bs="re")  +
                      s(plotID, bs="re") ,
                  family=nb, method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
summary(gam_cartae_re)
# A: plot and year

# A single common smoother plus group-level smoothers that have the same wiggliness (model GS)
gam_cartae_GS <- gam(sp_abund ~ s(DOY, k=5, bs="cc") +
                         s(collectDate, bs="re") + 
                         s(col_year_fac, bs="re") +
                         s(plotID, bs="re")  +
                         s(plot_trap, bs="re"),
                     family=nb, method="REML",
                     data=all7sp_dat %>% 
                         filter(para_sciname == "Carabus taedatus",
                                plotID != "NIWO_013"))
summary(gam_cartae_GS) 

gam_cartae <- gam(sp_abund ~ s(DOY, bs="cc", k=5) + te(DOY,precip_2weeks) + 
                      te(plot17aspect, plot17slope) +
                      te(LAI_1718avg, trap_CHM) + nlcdClass +
                        
                      s(collectDate, bs="re") + 
                      s(col_year_fac, bs="re") + 
                      s(plot_trap, bs="re")  +
                      s(plotID, bs="re") ,
                  family=nb, method="REML"
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"))
summary(gam_cartae)
jagam_cartae <- jagam(sp_abund ~ s(DOY, bs="cc", k=5) + te(DOY,precip_2weeks) + 
                      te(plot17aspect, plot17slope) +
                      te(LAI_1718avg, trap_CHM) + nlcdClass +
                        
                      s(collectDate, bs="re") + 
                      s(col_year_fac, bs="re") + 
                      s(plot_trap, bs="re")  +
                      s(plotID, bs="re") ,
                  family=poisson, 
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"),
                  file="jagam_cartae.txt")
summary(jagam_cartae)
jags_out <- jags()
gam_cartae_brm <- gam(sp_abund ~ s(DOY, bs="cc", k=5) + t2(DOY,precip_2weeks) + 
                      t2(plot17aspect, plot17slope) +
                      t2(LAI_1718avg, trap_CHM) + nlcdClass +
                        
                      s(collectDate, bs="re") + 
                      s(col_year_fac, bs="re") + 
                      s(plot_trap, bs="re")  +
                      s(plotID, bs="re") ,
                  family=poisson, 
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"),
                  file="jagam_cartae.txt")
summary(gam_cartae_brm)
brm(formula=bf(gam_cartae), 
    data=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"),
    family=negbinomial, 
    seed=2)


gam.check(gam_cartae) #check if k values are too low
summary(gam_cartae)
plot(gam_cartae)
draw(gam_cartae)
predict(gam_cartae, type = "response", newdata=all7sp_dat %>% 
                      filter(para_sciname == "Carabus taedatus",
                             plotID == "NIWO_013"))

# model evaluation: RMSE, mean absolute error, look at Bias (high or low?), or maybe it has good fit over certain time periods (better at winter months than summer, for example)
# see bird paper for model evaluation tips





# Cymindis unicolor
# What grouping variables explain the most variance? 
gam_cymuni_re <- gam(sp_abund ~ s(collectDate, bs="re") + 
                      s(collectDate, col_year_fac, bs="re") + 
                      s(plot_trap, bs="re")  +
                      s(plot_trap, plotID, bs="re") ,
                  family=nb, method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
summary(gam_cymuni_re)
# A: trap. Collectdate
gam_cymuni <- gam(sp_abund ~ s(DOY, bs="cc", k=3) + te(DOY, precip_2weeks) + te(plot17aspect,plot17slope) +
                      LAI_1718avg + trap_CHM + nlcdClass +
                      s(plotID, bs="re") + s(plot_trap, bs="re") + 
                      s(col_year_fac, bs="re") + s(collectDate, bs="re"), #grouping variables
                  family=nb(), method="REML",
                  data=all7sp_dat %>% 
                      filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"))
gam.check(gam_cymuni)
summary(gam_cymuni)
draw(gam_cymuni)
#

































# what's the difference between gam vs gamm fns?
# why smooth the grouping variables? ways to put grouping variables into a gam?
# compare a simulation of the fitted model to the data (predict)
# predict just predicts the means. then take random draws from the distribution
# does mgcv have a simulate fn?
# draw from the fitted mean, then from each grouping variable's normal dist, then from the final poisson dist
# plot simulated abund and actual abund vs each predictor variable
# Plot predictions of NIWO_0013





#for one species
# remove one plot or trap
# compare a few gams
#do model selection
# predict on removed trap/plot






# For Cymindis unicolor, we see higher abundance in the forest 
mod_cymuni_all <- gam(sp_abund ~ s(trap_LAI) +
                      trap_CHM +
                       as.numeric(as.character(col_year)) +
                       s(DOY,  bs = "cc", k=3) +
                       #s(elevation) +
                       #nlcdClass + #B/K says we want this to be a fixed effect
                       s(plotID, bs="re") +
                       s(plot_trap, bs="re") +
                       s(collectDate, bs="re"),
                   family=poisson(),
                   data=all7sp_dat %>% filter(para_sciname == "Cymindis unicolor"),
                   method="REML")
summary(mod_cymuni_all)
plot(mod_cymuni_all, pages=1)




# GLMs
mod_cartae_all <- glmer(sp_abund ~ trap_LAI +
                            trap_CHM +
                            col_year +
                            DOY^2 +
                            (1|plotID) +
                            (1|plot_trap) +
                            (1|collectDate),
                   family=poisson(),
                   data=all7sp_dat %>% filter(para_sciname == "Carabus taedatus"),
                   method="REML")


