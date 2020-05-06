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
    dplyr::select(elev, LAI_1718avg, plot_CHM, LAI_1718avg, precip_2weeks,
                  DOY, plot17aspect, plot17slope, GDD_cum, col_year, nlcdclass_num)
omcdiag(vars_mat, all7sp_dat$sp_abund)
imcdiag(vars_mat, all7sp_dat$sp_abund, method="VIF") #failed to dect multicollinearity when elevation was remoted

# remove elevation
vars_mat <- all7sp_dat %>%
    mutate(nlcdclass_num = as.numeric(nlcdClass)) %>%
    dplyr::select(LAI_1718avg, plot_CHM, LAI_1718avg, precip_2weeks, DOY,
                  plot17aspect, plot17slope, GDD_cum, col_year, nlcdclass_num)
omcdiag(vars_mat, all7sp_dat$sp_abund)
imcdiag(vars_mat, all7sp_dat$sp_abund, method="VIF") # when elev is removed, collinearity among other variables disappears.
# However, now plot_CHM, DOY, GDD_cum, nlcdclass_num have high VIF values

# Visualize correlation 
pairs(vars_mat, col="dodgerblue") 
# GDD_cum and DOY are nearly perfectly correlated. 

# Since GDD_cum is more biologically meaningful than DOY, drop DOY
vars_mat <- all7sp_dat %>%
    mutate(nlcdclass_num = as.numeric(nlcdClass)) %>%
    dplyr::select(LAI_1718avg, plot_CHM, LAI_1718avg, precip_2weeks, 
                  plot17aspect, plot17slope, GDD_cum, col_year, nlcdclass_num)
omcdiag(vars_mat, all7sp_dat$sp_abund)
imcdiag(vars_mat, all7sp_dat$sp_abund, method="VIF")
# nlcdclass_num and plot_CHM still have high VIF values. Try dropping plot_CHM, since it's similar to LAI and we can also model the other cont. variables by the categorical nlcdClass_num variable

vars_mat <- all7sp_dat %>%
    mutate(nlcdclass_num = as.numeric(nlcdClass)) %>%
    dplyr::select(LAI_1718avg, LAI_1718avg, precip_2weeks, 
                  plot17aspect, plot17slope, GDD_cum, col_year, nlcdclass_num)
omcdiag(vars_mat, all7sp_dat$sp_abund)
imcdiag(vars_mat, all7sp_dat$sp_abund, method="VIF")
# VIF values look good! 

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



# Model selection ---------------------------------------------------------

vp_mod_cartae <- gam(sp_abund ~ s(collectDate, bs="re") + 
                         s(col_year_fac, bs="re") + 
                         s(plot_trap, bs="re")  + s(plotID, bs="re") ,
                     family=nb, method="REML", data = all7sp_dat %>% 
                         filter(para_sciname == "Carabus taedatus", 
                                plotID!="NIWO_013"))
save(vp_mod_cartae, file="data_derived/vp_mod_cartae.Rdata")

vp_mod_cymuni <- gam(sp_abund ~ s(collectDate, bs="re") + 
                         s(col_year_fac, bs="re") + 
                         s(plot_trap, bs="re")  + s(plotID, bs="re") ,
                     family=nb, method="REML", data = all7sp_dat %>% 
                         filter(para_sciname == "Cymindis unicolor",
                                plotID!="NIWO_013"))
save(vp_mod_cymuni, file="data_derived/vp_mod_cymuni.Rdata")

for (spec in c("Carabus taedatus","Cymindis unicolor")) {
    
    #Null model
    null_mod <- gam(sp_abund ~ 1, family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    
    # Variance partition model
    vp_mod <- gam(sp_abund ~ s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                      s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    
    # Negative binomial models
    # Just nlcdClass
    nb2_mod <- gam(sp_abund ~ nlcdClass +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Just slope/aspect interaction
    nb3_mod <- gam(sp_abund ~ te(plot17aspect, plot17slope) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Just LAI
    nb4_mod <- gam(sp_abund ~ s(LAI_1718avg) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Just CHM
    nb5_mod <- gam(sp_abund ~ s(trap_CHM) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Just precip
    nb6_mod <- gam(sp_abund ~ s(precip_2weeks) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Just GDD
    nb7_mod <- gam(sp_abund ~ s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Just DOY
    nb8_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    
    # Add nlcdClass
    nb9_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add slope/aspect interaction
    nb10_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add precip
    nb11_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add precp x DOY interaction
    nb12_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, bs=c("cc","ts")) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add LAI
    nb13_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, bs=c("cc","ts")) + s(LAI_1718avg) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add CHM - shouldn't perform well since ther are highly correlated variables
    nb14_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, bs=c("cc","ts")) + s(LAI_1718avg) + s(plot_CHM, k=3) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add GDD - shouldn't perform well since ther are highly correlated variables
    nb15_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, bs=c("cc","ts")) + s(LAI_1718avg) + s(plot_CHM, k=3) + s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Adjust for highly correlated variables: canopy height, DOY
    nb16_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, bs=c("cc","ts")) + s(LAI_1718avg) + s(plot_CHM, k=3) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Adjust for highly correlated variables: canopy height, GDD
    nb17_mod <- gam(sp_abund ~ te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks) + s(LAI_1718avg) + s(plot_CHM, k=3) + s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Adjust for highly correlated variables: nlcdClass, DOY
    nb18_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, bs=c("cc","ts")) + s(LAI_1718avg) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Adjust for highly correlated variables: nlcdClass, GDD
    nb19_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, bs=c("cc","ts")) + s(LAI_1718avg) +s( GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Replace LAI and canopy height with interaction
    nb20_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks) + te(LAI_1718avg,plot_CHM) + s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Try variables that stepGAM said were important
    nb21_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks) + s(LAI_1718avg) + s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Try variables that stepGAM said were important - no GDD_cum
    nb22_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        te(DOY,precip_2weeks) + s(LAI_1718avg) + 
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add interaction between precip and GDD_cum
    nb23_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope) + s(precip_2weeks) +
                        s(LAI_1718avg) + s(GDD_cum) + te(precip_2weeks, GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add interaction between precip and GDD_cum
    nb24_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope) + 
                        te(precip_2weeks, GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Try nlcdClass, slope/aspect, CHM, DOY/precip
    nb25_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope) + 
                        s(plot_CHM, k=3) + te(DOY, precip_2weeks, bs=c("cc","ts")) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    
    # Now, allow intercept to vary by nlcdClass by adding by=nlcdClass to each plot- or trap-level continuous variable
    # Adjust for highly correlated variables: canopy height, DOY
    nb26_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + te(plot17aspect, plot17slope, by=nlcdClass) + 
                        s(precip_2weeks, by=nlcdClass) +
                        te(DOY,precip_2weeks, bs=c("cc","ts"), by=nlcdClass) + s(LAI_1718avg, by=nlcdClass) + s(plot_CHM, by=nlcdClass, k=3) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Adjust for highly correlated variables: canopy height, GDD
    nb27_mod <- gam(sp_abund ~ te(plot17aspect, plot17slope, by=nlcdClass) + s(precip_2weeks, by=nlcdClass) +
                        te(DOY,precip_2weeks) + s(LAI_1718avg, by=nlcdClass) + s(plot_CHM, by=nlcdClass, k=3) + s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Adjust for highly correlated variables: nlcdClass, DOY
    nb28_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + s(precip_2weeks) +
                        te(DOY,precip_2weeks, by=nlcdClass, bs=c("cc","ts")) + s(LAI_1718avg, by=nlcdClass) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Adjust for highly correlated variables: nlcdClass, GDD
    nb29_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + s(precip_2weeks, by=nlcdClass) +
                        te(DOY,precip_2weeks, by=nlcdClass, bs=c("cc","ts")) + s(LAI_1718avg, by=nlcdClass) +s( GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Replace LAI and canopy height with interaction
    nb30_mod <- gam(sp_abund ~ s(DOY, bs="cc",k=5) + nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + s(precip_2weeks, by=nlcdClass) +
                        te(DOY,precip_2weeks, by=nlcdClass, bs=c("cc","ts")) + te(LAI_1718avg,plot_CHM) + s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Try variables that stepGAM said were important
    nb31_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + s(precip_2weeks, by=nlcdClass) +
                        te(DOY,precip_2weeks, by=nlcdClass, bs=c("cc","ts")) + s(LAI_1718avg, by=nlcdClass) + s(GDD_cum) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Try variables that stepGAM said were important - no GDD_cum
    nb32_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + s(precip_2weeks, by=nlcdClass) +
                        te(DOY,precip_2weeks, by=nlcdClass, bs=c("cc","ts")) + s(LAI_1718avg, by=nlcdClass) + 
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add interaction between precip and GDD_cum
    nb33_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + s(precip_2weeks, by=nlcdClass) +
                        s(LAI_1718avg, by=nlcdClass) + s(GDD_cum) + te(precip_2weeks, GDD_cum, by=nlcdClass) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Add interaction between precip and GDD_cum
    nb34_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + 
                        te(precip_2weeks, GDD_cum, by=nlcdClass) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    # Try nlcdClass, slope/aspect, CHM, DOY/precip
    nb35_mod <- gam(sp_abund ~ nlcdClass + te(plot17aspect, plot17slope, by=nlcdClass) + 
                        s(plot_CHM, by=nlcdClass, k=3)+ te(precip_2weeks, DOY, by=nlcdClass, bs=c("cc","ts")) +
                       s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                    s(plot_trap, bs="re")  + s(plotID, bs="re") , family=nb, method="REML",
                     data=all7sp_dat %>% filter(para_sciname == spec, plotID != "NIWO_013"))
    
    # stop there, but then go back to add by=nlcdClass
    

    
    # Zero-inflated Poisson models
    
    
    # Save models
    pattern <- ls(pattern="_mod") # get names of models
    models <- do.call("list", mget(pattern))
    if (spec == "Carabus taedatus") {
        save(models, file="data_derived/model_sel_cartae.Rdata", overwrite=T)
    } else {
        save(models, file="data_derived/model_sel_cymuni.Rdata", overwrite=T)
    }
}





all7sp_dat <- read.csv("data_derived/model_df_by_species_in_sample.csv") %>%
    mutate(sc_DOY = scale(DOY, center = TRUE, scale = TRUE),
           col_year_fac = as.factor(col_year))

# Load in models 
load("data_derived/model_sel_cartae.Rdata")
models_cartae <- models
for (i in 1:length(models_cartae)) {
    names(models_cartae)[i] <- paste0(names(models_cartae)[i],"_cartae",sep="")
}
load("data_derived/model_sel_cymuni.Rdata")
models_cymuni <- models
for (i in 1:length(models_cymuni)) {
    names(models_cymuni)[i] <- paste0(names(models_cymuni)[i],"_cymuni",sep="")
}
rm(models)
list2env(models_cartae, .GlobalEnv)
list2env(models_cymuni, .GlobalEnv)
# do.call(rbind, lapply(models_cartae, glance)) %>% 
#     data.frame() %>%
#     mutate(model = pattern) %>% 
#     select(model, everything(.)) %>% 
#     select(-c(df.residual, BIC, logLik)) %>%
#     arrange(-AIC)
# nb4_mod_cartae$deviance
# nb4_mod_cartae$null.deviance #not sure what this is
# nb4_mod_cartae$aic # different from the aic reported above
# nb4_mod_cartae$sp #lists random effects
# nb4_mod_cartae$rank #not sure what this is
# nb4_mod_cartae$gcv.ubre #prints REML
# nb4_mod_cartae$formula #lists formula or $call
# nb4_mod_cartae$df.residual#good metric


# We compared a ton of models above and looked to see which predictors were most useful for each species. Now, let's take those most imporant predictors and create model iterations
    # Try nlcdClass, slope/aspect, CHM, DOY/precip
dat_cartae <- all7sp_dat %>% filter(para_sciname == "Carabus taedatus", plotID != "NIWO_013")

glob_nb_mod_cartae <- gam(sp_abund ~ nlcdClass + s(DOY,bs="cc",k=3) + s(precip_2weeks, bs="ts") +
                              s(precip_2weeks, by=nlcdClass, bs="ts") +
                              s(plot_CHM, bs="ts",k=5) + s(plot_CHM, by=nlcdClass, bs="ts",k=5) +
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)

# Only with 'by-nlcdClass'
nb1_mod_cartae <- gam(sp_abund ~ nlcdClass + s(DOY,bs="cc",k=3) + 
                          s(precip_2weeks, by=nlcdClass, bs="ts") +
                          s(plot_CHM, by=nlcdClass, bs="ts",k=5) +
                          te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                          te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
# Without 'by-nlcdClass'
nb2_mod_cartae <- gam(sp_abund ~ nlcdClass + s(DOY,bs="cc",k=3) + 
                          s(precip_2weeks, bs="ts") +
                          s(plot_CHM, bs="ts",k=5) + s(plot_CHM, by=nlcdClass, bs="ts",k=5) +
                          te(DOY, precip_2weeks, bs=c("cc","ts")) +
                          te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)

# Does te vs s matter for slope/aspect?
nb3_mod_cartae <- gam(sp_abund ~ nlcdClass + s(DOY,bs="cc",k=3) + 
                          s(precip_2weeks, bs="ts",k=3) +
                          s(plot_CHM, bs="ts",k=3) + s(plot_CHM, by=nlcdClass, bs="ts",k=3) +
                          te(DOY, precip_2weeks, bs=c("cc","ts"),k=3) +
                          s(plot17aspect, plot17slope, bs=c("ts","ts"),k=3) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)

# Look at summary output and then remove some
nb4_mod_cartae <- gam(sp_abund ~ nlcdClass + s(DOY,bs="cc",k=3) + 
                              s(precip_2weeks, by=nlcdClass, bs="ts") +
                              s(plot_CHM, bs="ts",k=5) + s(plot_CHM, by=nlcdClass, bs="ts",k=5) +
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
nb5_mod_cartae <- gam(sp_abund ~ nlcdClass + s(DOY,bs="cc",k=3) + 
                              s(precip_2weeks, by=nlcdClass, bs="ts") +
                              s(plot_CHM, bs="ts",k=5) + 
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
nb6_mod_cartae <- gam(sp_abund ~ nlcdClass + 
                              s(precip_2weeks, by=nlcdClass, bs="ts") +
                              s(plot_CHM, bs="ts",k=5) + 
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
nb7_mod_cartae <- gam(sp_abund ~ nlcdClass + 
                              s(plot_CHM, bs="ts",k=5) + 
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
nb8_mod_cartae <- gam(sp_abund ~ nlcdClass + 
                              s(plot_CHM, bs="ts",k=5) + 
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
# nb5 without nlcdClass
nb9_mod_cartae <- gam(sp_abund ~ s(DOY,bs="cc",k=3) + 
                              s(precip_2weeks, by=nlcdClass, bs="ts") +
                              s(plot_CHM, bs="ts",k=5) + 
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)

# All of the models are within 2 AIC values, so let me start to try to differentiate
nb10_mod_cartae <- gam(sp_abund ~ nlcdClass + s(DOY,bs="cc",k=3) + 
                              s(precip_2weeks, by=nlcdClass, bs="ts") +
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
nb11_mod_cartae <- gam(sp_abund ~ nlcdClass + s(GDD_cum) + 
                              s(precip_2weeks, by=nlcdClass, bs="ts") +
                              te(DOY, precip_2weeks, bs=c("cc","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)
nb12_mod_cartae <- gam(sp_abund ~ nlcdClass + s(GDD_cum) + 
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts")) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                          s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                          s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                      family=nb, method="REML", data=dat_cartae)

model.sel(glob_nb_mod_cartae, nb1_mod_cartae, nb2_mod_cartae, nb3_mod_cartae, nb4_mod_cartae, nb5_mod_cartae, nb6_mod_cartae, nb7_mod_cartae, nb8_mod_cartae, nb9_mod_cartae,nb10_mod_cartae, nb11_mod_cartae, nb12_mod_cartae, rank=AIC)

# So many of the models are equivalent by AIC, so arbitrarily choose one
#save(nb5_mod_cartae, file="data_derived/final_mod_cartae.Rdata")

#AIS another time, try dredge for exhaustive model comparison. Took WAY too long to process to be usable for our presentation (right now, it's been 6hrs): dredge_cartae <- lapply(dredge(glob_nb_mod_cartae, evaluate = FALSE), eval)

dat_cymuni <- all7sp_dat %>% filter(para_sciname == "Cymindis unicolor", plotID != "NIWO_013")
glob_nb_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,bs="ts") + s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,bs="ts") + s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks,bs="ts", k=5) + s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, bs=c("cc","ts"), k=4) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, bs=c("ts","ts"), k=3) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)

glob_nonlcd_nb_mod_cymuni <- gam(sp_abund ~ 
                              s(GDD_cum,bs="ts") + s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,bs="ts") + s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks,bs="ts", k=5) + s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, bs=c("cc","ts"), k=4) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, bs=c("ts","ts"), k=3) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)

# Only 'by nlcdclass' terms
nb1_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)

# No 'by nlcdclass' terms
nb2_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,bs="ts") +
                              s(LAI_1718avg,bs="ts") + 
                              s(precip_2weeks,bs="ts", k=5) + 
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, bs=c("ts","ts"), k=3) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)

# Half glob - top half with nlcdclass
halfglob1_nb_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, bs=c("ts","ts"), k=3) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
# Half glob - bottom half with nlcdclass
halfglob2_nb_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,bs="ts") + 
                              s(LAI_1718avg,bs="ts") + 
                              s(precip_2weeks,bs="ts", k=5) + 
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
model.sel(halfglob1_nb_mod_cymuni, halfglob2_nb_mod_cymuni, rank=AIC)

# mod1 had the lowest AIC. Check out various iterations of it
nb1.1_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb1.2_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb1.3_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb1.4_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb1.5_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb1.6_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb1.7_mod_cymuni <- gam(sp_abund ~ 
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)

# Try a few other variations
nb3_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks,bs="ts", k=5) + s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, bs=c("cc","ts"), k=4) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, bs=c("ts","ts"), k=3) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb4_mod_cymuni <- gam(sp_abund ~ 
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb5_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,by=nlcdClass,bs="ts") +
                              s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, GDD_cum, by=nlcdClass, bs=c("cc","ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb6_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(LAI_1718avg,bs="ts") + s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)

# removing predictors after looking at summary and plots of global model
nb7_mod_cymuni <- gam(sp_abund ~ #nlcdClass +
                              s(GDD_cum,bs="ts") + 
                              LAI_1718avg +
                              precip_2weeks +
                              #s(LAI_1718avg,bs="ts") + s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              # summary says LAI and precip were useful, but it jsut doesn't look right
                              #s(precip_2weeks,bs="ts", k=5) + s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, bs=c("cc","ts"), k=4) +
                              te(DOY, precip_2weeks, by=nlcdClass, bs=c("cc","ts"), k=4) +
                              te(GDD_cum, precip_2weeks, bs=c("ts","ts"), k=3) +
                              te(GDD_cum, precip_2weeks, by=nlcdClass, bs=c("ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb8_mod_cymuni <- gam(sp_abund ~ #nlcdClass +
                              s(GDD_cum,bs="ts") + 
                              LAI_1718avg +
                              precip_2weeks +
                              #s(LAI_1718avg,bs="ts") + s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              # summary says LAI and precip were useful, but it jsut doesn't look right
                              #s(precip_2weeks,bs="ts", k=5) + s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, GDD_cum, bs=c("cc","ts","ts"), k=4) +
                              te(DOY, precip_2weeks, GDD_cum, by=nlcdClass, bs=c("cc","ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)
nb9_mod_cymuni <- gam(sp_abund ~ nlcdClass +
                              s(GDD_cum,bs="ts") + 
                              LAI_1718avg +
                              precip_2weeks +
                              #s(LAI_1718avg,bs="ts") + s(LAI_1718avg,by=nlcdClass,bs="ts") +
                              # summary says LAI and precip were useful, but it jsut doesn't look right
                              #s(precip_2weeks,bs="ts", k=5) + s(precip_2weeks, by=nlcdClass,bs="ts", k=5) +
                              te(plot17aspect, plot17slope, bs=c("ts","ts")) +
                              te(plot17aspect, plot17slope, by=nlcdClass, bs=c("ts","ts")) +
                              te(DOY, precip_2weeks, GDD_cum, bs=c("cc","ts","ts"), k=4) +
                              te(DOY, precip_2weeks, GDD_cum, by=nlcdClass, bs=c("cc","ts","ts"), k=4) +
                              s(collectDate, bs="re") + s(col_year_fac, bs="re") + 
                              s(plot_trap, bs="re")  + s(plotID, bs="re"), 
                          family=nb, method="REML", data=dat_cymuni)


model.sel(glob_nb_mod_cymuni, glob_nonlcd_nb_mod_cymuni, nb1_mod_cymuni, nb2_mod_cymuni, nb1_mod_cymuni, nb1.1_mod_cymuni, nb1.2_mod_cymuni, nb1.3_mod_cymuni, nb1.4_mod_cymuni, nb1.5_mod_cymuni, nb1.6_mod_cymuni, nb1.7_mod_cymuni, halfglob1_nb_mod_cymuni, halfglob2_nb_mod_cymuni, nb3_mod_cymuni, nb4_mod_cymuni, nb5_mod_cymuni, nb6_mod_cymuni, nb7_mod_cymuni, nb8_mod_cymuni, nb9_mod_cymuni, rank=AIC) 

final_mods_cymuni <- list(glob_nb_mod_cymuni, glob_nonlcd_nb_mod_cymuni, nb1.7_mod_cymuni, nb1.4_mod_cymuni)
#save(final_mods_cymuni, file="data_derived/final_mods_cymuni.Rdata")
    

# Final models - convert te to t2 -----------------------------------------

load("data_derived/final_mod_cartae.Rdata")
load("data_derived/final_mod_cymuni.Rdata")
cartae_mod <- nb5_mod_cartae
cymuni_mod <- final_mod_cymuni
rm(nb5_mod_cartae, final_mod_cymuni)

cartae_mod$formula
cartae_mod_t2 <- gam(sp_abund ~ nlcdClass + s(DOY, bs = "cc", k = 3) + s(precip_2weeks, 
    by = nlcdClass, bs = "ts") + s(plot_CHM, bs = "ts", k = 5) + 
    t2(DOY, precip_2weeks, bs = c("cc", "ts")) + t2(DOY, precip_2weeks, 
    by = nlcdClass, bs = c("cc", "ts")) + t2(plot17aspect, plot17slope, 
    bs = c("ts", "ts")) + t2(plot17aspect, plot17slope, by = nlcdClass, 
    bs = c("ts", "ts")) + s(collectDate, bs = "re") + s(col_year_fac, 
    bs = "re") + s(plot_trap, bs = "re") + s(plotID, bs = "re"), data =  all7sp_dat %>% filter(para_sciname == "Carabus taedatus", plotID != "NIWO_013"), family = nb, method="REML")
save(cartae_mod_t2, file="data_derived/final_mod_cartae_t2.Rdata")

cymuni_mod$formula
cymuni_mod_t2 <- gam(sp_abund ~ nlcdClass + s(GDD_cum, bs = "ts") + s(GDD_cum, by = nlcdClass, bs = "ts") + s(LAI_1718avg, bs = "ts") + s(LAI_1718avg, by = nlcdClass, 
    bs = "ts") + s(precip_2weeks, bs = "ts", k = 5) + s(precip_2weeks, 
    by = nlcdClass, bs = "ts", k = 5) + t2(plot17aspect, plot17slope, 
    bs = c("ts", "ts")) + t2(plot17aspect, plot17slope, by = nlcdClass, 
    bs = c("ts", "ts")) + t2(DOY, precip_2weeks, bs = c("cc", 
    "ts"), k = 4) + t2(DOY, precip_2weeks, by = nlcdClass, bs = c("cc", 
    "ts"), k = 4) + t2(GDD_cum, precip_2weeks, bs = c("ts", "ts"), 
    k = 3) + t2(GDD_cum, precip_2weeks, by = nlcdClass, bs = c("ts", 
    "ts"), k = 4) + s(collectDate, bs = "re") + s(col_year_fac, 
    bs = "re") + s(plot_trap, bs = "re") + s(plotID, bs = "re"), data =  all7sp_dat %>% filter(para_sciname == "Cymindis unicolor", plotID != "NIWO_013"), family = nb, method="REML")
save(cymuni_mod_t2, file="data_derived/final_mod_cymuni_t2.Rdata")

# Do they differ in AIC?
model.sel(cartae_mod, cartae_mod_t2)
model.sel(cymuni_mod, cymuni_mod_t2)


# AIS, when I return to this, I should set up an nxn table where n is the number of possible predictors (including interactions) and then I need to iterate through all of the unique combinations of these



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


# smooth everythin
# bs = ts (shrinkage) vs tp (default) - what does wynne do?
# include shrinkage. If a term is using effectively 0 degrees of freedom it is having no effect on the fit/predictions at all. For the non-significant terms that have positive EDFs, by keeping them in you are effectively stating that these covariates have a small but non-zero effect. If you remove these terms as you suggest, you are saying explicitly that the effect is zero.

# Try Melissa's stepGAM fn
allPred <- c("s(DOY, bs='cc', k=3)", "s(GDD_cum, bs='ts')", "te(DOY, precip_2weeks)", 
             "te(plot17aspect,plot17slope)", "s(LAI_1718avg, bs='ts')", "s(trap_CHM, bs='ts')", "nlcdClass",
                      "s(plotID, bs='re')", "s(plot_trap, bs='re')", 
                      "s(col_year_fac, bs='re')", "s(collectDate, bs='re')")
allAIC_cymuni <- stepGAM(dat = all7sp_dat %>%
                             filter(para_sciname == "Cymindis unicolor",
                             plotID != "NIWO_013"),
                         predictors = allPred, 
                         response = "sp_abund", 
                         family = nb,
                         ignore.combos=NA)
save(allAIC_cymuni, file="data_derived/stepGAM_output_cymuni.Rdata")
allAIC_cartae <- stepGAM(dat = all7sp_dat %>%
                             filter(para_sciname == "Carabus taedatus",
                             plotID != "NIWO_013"),
                         predictors = allPred, 
                         response = "sp_abund", 
                         family = nb,
                         ignore.combos=NA)
save(allAIC_cartae, file="data_derived/stepGAM_output_cartae.Rdata")

