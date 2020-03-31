# We try different models of predicting carabid abundance and composition

library(lme4)      #max lik multilevel: lmer(), glmer() etc
library(arm)       #for se.ranef()
library(ggplot2)
library(gridExtra) #arranging multiple plots
library(dplyr)
library(mgcv) #for gams

all7sp_dat <- read.csv("data_derived/model_df_by_species_in_sample.csv")

# Pick one species for now. 
all7sp_dat %>% group_by(para_sciname) %>%
    summarize(total = sum(sp_abund)) %>%
    arrange(-total)
# Calathus advena is most abundant
caladv_dat <- all7sp_dat %>% 
    filter(para_sciname == "Calathus advena") %>%
    arrange(collectDate) %>%
    mutate_at(c("siteID", "plotID", "trapID", "sampleID", "col_year", "col_month", "col_day", "nlcdClass", "soilOrder", "para_sciname"), funs(factor(.)))
    

# Try modeling Calathus advena
# Following Brett's multilevel tutorials

# EDA -------------------------

# How many samples with 0 abundance?
sum(caladv_dat$sp_abund == 0) 
# 813 out of 960
ggplot(data=caladv_dat) +
    geom_histogram(mapping = aes(x=sp_abund),bins=30)
# Heavy 0 count data

# Try log-abundance
ggplot(data=caladv_dat) +
    geom_histogram(mapping = aes(x=log(sp_abund),bins=30))
# looks pretty similar

ggplot(data=caladv_dat) +
    geom_histogram(mapping = aes(x=sp_abund, fill=nlcdClass), position="identity",
                   bins=20,alpha=0.5)
# This species is found only in forest plots, not on the tundra

# Plot 4 was swapped for plot 13 in 2018. Compare plot abundance in first 3 years to see if plot 4 was an outlier
all7sp_dat %>%
    ggplot() +
    geom_point(mapping = aes(x=dayofyear,col=col_year,y=sp_abund),shape=1,alpha=0.75) +
    facet_wrap(facets = ~ plotID) 
# Looking at 7 most abundant species, it seems that plot 4 didn't have have large beetle yield.... maybe that's why they moved it?
caladv_dat %>%
    group_by(plotID) %>%
    summarize(abund = sum(sp_abund)) %>%
    arrange(-abund)
#Calathus advena is present in only 7 of the 11 plots
caladv_dat %>%
    ggplot() +
    geom_point(mapping = aes(x=dayofyear, colour=as.numeric(as.character(col_year)), y=sp_abund),shape=1) +
    facet_wrap(facets = ~ plotID) 
# Consider dropping plot 4 from the data?

# Density
caladv_dat %>%
    ggplot() +
    geom_density(mapping = aes(x=sp_abund, col=col_month))
# All years had proportionally high 0 counts

# Through day of year
caladv_dat %>%
    ggplot(aes(x=dayofyear, y=sp_abund)) +
    geom_bar(stat="identity") +
    facet_grid(col_year ~ .) 
# Not really regular abundance patterns throughout the year.... not sure how to summarize



### Multilevel analysis --------------
# including the many 0-count samples

# Complete pooling
poolmean <- mean(caladv_dat$sp_abund)
poolmean
cp_pred_df <- data.frame(poolmean) #df for use with ggplot

# No pooling
# Does this make sense to do since abundance is not a normally distributed outcome varable
sp_abund_mean_var <- 
    caladv_dat %>%
    group_by(plotID) %>%
    summarize(sample_size=n(),plt_mn=mean(sp_abund),plt_sd=sd(sp_abund)) %>%
    mutate(plt_se=plt_sd/sqrt(sample_size)) %>%
    mutate(sample_size_jit=jitter(sample_size)) #jitter added for plotting
print(sp_abund_mean_var,n=Inf) #n=Inf to print all rows

ggplot(data=sp_abund_mean_var) +
    geom_hline(mapping=aes(yintercept=poolmean),data=cp_pred_df,col="blue") +
    geom_point(mapping=aes(x=sample_size_jit,y=plt_mn)) +
    geom_text(aes(x=sample_size_jit,y=plt_mn,label=plotID), hjust=-.1, vjust=0) +
    geom_linerange(mapping=aes(x=sample_size_jit,ymin=plt_mn-plt_se,ymax=plt_mn+plt_se)) +
    #scale_x_continuous(trans="log",breaks=c(1,3,10,30,100)) +
    labs(x="Sample size in plot j",y="mean abundance of Calathus advena in plot j",
         title="No pooling: separate analyses by plot")

# Partial pooling - plot as fixed effect
# Try GLM using poisson distribution for abundance
npfit_pois <- glm( sp_abund ~ -1 + plotID, family=poisson, data=caladv_dat )
plot(npfit_pois,1:5,ask=FALSE)
# This looks cray

#GAMs
caladv_dat_gam <- caladv_dat %>%
    mutate(plot_trap = as.factor(paste(plotID, trapID, sep="")),
           occ = ifelse(sp_abund>0,1,0),
           scaled_elev = scale(elevation, center = TRUE, scale = TRUE),
           scaled_dayofyear = scale(dayofyear, center = TRUE, scale = TRUE))

modre <- gam(sp_abund ~ s(nlcdClass, bs="re"),
             s(plot_trap, bs="re"),
             data=caladv_dat_gam, method="REML")

mod1 <- gam(sp_abund ~ s(elevation),
            s(nlcdClass, bs="re"),
            s(plot_trap, bs="re"),
           data=caladv_dat_gam, method="REML")

mod2d_occ <- gam(occ ~ s(decimalLatitude, decimalLongitude),
            #s(elevation),
            #s(nlcdClass,bs="re"),
            s(plotID, bs="re"),
            s(plot_trap, bs="re"),
            family=poisson(),
            data=caladv_dat_gam, method="REML")
mod2d_abund <- gam(sp_abund ~ s(as.numeric(as.character(col_year)), k=3) +
                       s(dayofyear,  bs = "cc", k=3) +
                       s(elevation) +
                       s(nlcdClass,bs="re") +
                       s(plotID, bs="re") +
                       s(plot_trap, bs="re"),
                   family=poisson(),
                   data=caladv_dat_gam, method="REML")
plot(mod2d_abund, pages=1, all.terms=TRUE, scale=0)
summary(mod2d_abund)
#try gam with year term with other species - this one is speciose always, so may not fluxuate as much between 
# random effects are pretty much priors
#jSDMs are an attempt to incorporate biotic interactions, but this is debated
# find overwintering behavior for each species


mod2d_elev <- gam(sp_abund ~ s(elevation),
            s(nlcdClass,bs="re"),
            s(plotID, bs="re"),
            s(plot_trap, bs="re"),
            family=poisson(),
            data=caladv_dat_gam, method="REML")
summary(mod1)
plot(mod2d, scheme=2)




caladv_dat <- caladv_dat %>%
    mutate(scaled_elev = scale(elevation, center = TRUE, scale = TRUE),
           scaled_dayofyear = scale(dayofyear, center = TRUE, scale = TRUE))
var_part <-glmer( sp_abund ~ scaled_elev + scaled_dayofyear +
                      (1|plotID/trapID) + 
                      (1|nlcdClass), 
                  family=poisson, data=caladv_dat )
plot(var_part)
summary(var_part)
plot(var_part,ask=FALSE)
library(DHARMa)
plot(caladv_dat$sp_abund, predict(var_part, type="response"))
# nlcdcover tricky bc there are no beetles in herbaceous
# try different distribution (not poisson)
# could be temporal
# create a cyclic spline for dayofyear (mgcv, "cc" for cyclic, "re" for random effect) - with gamms, recommended to create trap_plotID variables to capture nested nature
var_part_spl 

# Try GLM using quasipoisson distribution for abundance
npfit_qpois <- glm( sp_abund ~ -1 + plotID, family=quasipoisson, data=caladv_dat )
plot(npfit_qpois,1:5,ask=FALSE)
# Still cray, but slightly better

# Try GLM using neg-bin distribution for abundance
npfit_nb <- glm.nb(sp_abund ~ -1 + plotID, data=caladv_dat)
plot(npfit_nb,1:5,ask=FALSE)
# Still cray, similar to quasipoisson

npfit_qpois_ln <- glm( log(sp_abund) ~ -1 + plotID, family=quasipoisson, data=caladv_dat )
plot(npfit_qpois_ln,1:5,ask=FALSE)

var_part <-glmer( sp_abund ~ 1 + (1|plotID/trapID), family=poisson, data=caladv_dat )
plot(var_part)
summary(var_part)







var_part <-glmer( sp_abund ~ 1 + (1|plotID/trapID), family=poisson, data=caladv_dat )
plot(var_part)
summary(var_part)

# Let's go with quasipoisson dist.
# Now plot the residuals
r <- residuals(var_part)
#r <- residuals(npfit_qpois_ln)
x <- seq(min(r),max(r),length.out=100)
y <- dnorm(x,mean(r),sd(r))
res_df <- data.frame(residuals=r)
norm_df <- data.frame(x=x,y=y)
rm(r,x,y)
ggplot() +
    geom_histogram(mapping=aes(x=residuals,y=stat(density)),data=res_df,bins=50) +
    geom_line(mapping=aes(x=x,y=y),col="red",data=norm_df)
# AIS what to make of this residual distribution?? many small negative residuals and few large postitive residuals. why would the residuals be lognormal/

np_pred_df <- data.frame(sp_abund_mean_var$plotID, coef(summary(var_part)),sp_abund_mean_var$sample_size_jit)
names(np_pred_df) <- c("plotID", "plt_mn","plt_se","sample_size_jit")
gh12.1a <- 
    ggplot(data=np_pred_df) +
    geom_hline(mapping=aes(yintercept=poolmean),data=cp_pred_df,col="blue") +
    geom_point(mapping=aes(x=sample_size_jit,y=plt_mn)) +
    geom_text(aes(x=sample_size_jit,y=plt_mn,label=plotID), hjust=-.1, vjust=0) +
    geom_linerange(mapping=aes(x=sample_size_jit,ymin=plt_mn-plt_se,ymax=plt_mn+plt_se)) +
    labs(x="Sample size in plot j",y="mean abundance of Calathus advena in plot j",
         title="No pooling: estimates from linear model fit")
gh12.1a
# That doesn't look right at all, but at least the plots with higher sample size have smaller SEs

# Partial pooling - plot as random effect
# Try without and with logging abundance
ppfit <- glmer( sp_abund ~ 1 + (1|plotID), REML=FALSE, family=poisson, data=caladv_dat )
plot(ppfit)
ppfit <- glmer( log(sp_abund) ~ 1 + (1|plotID), REML=FALSE, family=poisson, data=caladv_dat )
plot(ppfit)
# Logging abundance has fewer outliers, let's go with that

summary(ppfit)
# Hmm error. Try removing plot 4 rows, only 2 rows
subset(caladv_dat, plotID == "NIWO_004")
caladv_dat_wo4 <- subset(caladv_dat, plotID != "NIWO_004")

ppfit <- glmer( log(sp_abund) ~ 1 + (1|plotID), REML=FALSE, family=poisson, data=caladv_dat_wo4 )
plot(ppfit)

var_part <-glmer( sp_abund ~ 1 + (1|plotID/trapID), family=poisson, data=caladv_dat )
plot(var_part)
summary(var_part)


# prediction intervals
# cross-validation analysis
# see time series tutorial
# see timetable for next TO DO
# next time, see complete df







# Poisson, quasipoisson, or neg-bin?

# Harris' baseline average model
# Investivating tthe scales of variance in the data (aka variance components model). In Dietze book, he uses the average model as an opportunity of exploring the scales of variance in the data. 
# Between trap variation and between plot variation.
# then you have multiple samples per trap through time
# Where is most of the variation? Is it between plots or between traps?
# This average model with a random effects structure is farily sophisticated at different scales.

# Harris' baseline naive model
