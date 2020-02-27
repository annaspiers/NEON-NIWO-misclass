# We try different models of predicting carabid abundance and composition

library(arm)

dat <- lmer

# Try modeling

fit_abund <- glmer(formula = sp_abund ~ nlcdClass + para_sciname + (1|plotID) , 
    family = poisson(link="log"),
      data = dat)
display(fit_abund)