---
title: "Forecasting Carabid Beetle Abundance at Niwot Ridge"
author: "Anna Spiers, Christa Torrens, Grant Vagle"
date: "5/07/2020"
output:
  html_document:
    df_print: paged
    keep_md: true
---



# Project overview and background

## Project goals - big picture

* To use NEON data, including AOP lidar data, to accurately forecast the abundance of carabid beetles across the landscape, using training data from Niwot Ridge. 

* To be able to forecast changes to carabid beetle distribution under changing climate conditions 

* To use our forecast model to determine "how much data is enough?" for accurate forcasting at different timesteps
    + Assume the predicted data are real - run model 25 years into future
    + Test model accuracy with 1, 3, 5, 10, 15, 20 years of “real” data


## Project goals - this semester

* To create and validate a working model to accurately predict the abundance of 1-2 carabid spp. at NEON's Niwot carabid traps 
    + *either* a temporal prediction (2019 season)
    + *or* spatial prediction (one of the plots, 2018 season).


## Why Carabids?

Ground beetles (Carabidae) are an excellent sentinel taxon, and can be used to understand the consequences of global climate change in high mountain areas. 

* They are abundant, found across many biomes, and are relatively easy to identify

* They are sensitive to habitat change, thus are predicted to serve as useful bioindicators of environmental and land use change 
    + Numerous studies have shown that carabid beetle richness and/or assemblage composition change along altitudinal gradients, with vegetation type, and/ or land use type 
    
* Their presence/ absence impacts other ecosystem components: 
    + They are an important component of terrestrial food webs and can influence trophic structure 
(Hoekman et al. 2017, Hiramatsu and Usio 2018)

Available data: The NEON data suite includes carabid abundance since 2015, as well as many of our target predictor variables. 


## Research questions (initial)

1. Which carabid species are 'best' to focus on for our semester goals?

2. Initial exploration of variables (EDA)
+ Based on the literature, what are likely to be important variables driving carabid presence/ absence? 
    + Note: We initially focused on spatial predictions, but found temporal variance was a stronger influence on carabid abundance

3. What are the key predictor variables for each species?  (Model building)
    
4. Which model type will best predict species abundance? (Model testing)


## Q1: Carabid species selection

First we determined that there were seven abundant and well-identified species in NEON's Niwot data.  

<img src="images/Carabid ID_7circled.png" width="600" />

(the next section will go into details on sampling design and analysis)

Then we looked at abundance by sampling plot and collection date

<img src="output/species_through_years.png" width="3283" />

Looking at the abundance data, we chose two of those seven: *Carabus taedatus* and *Cymindis unicolor*.

<img src="output/abund_DateClass_twospp.png" width="1920" />

<!-- find code to create the distribution by nlcd class and spp image --> 

Both species were consistently present in both of Niwot's landscape classes: tundra and evergreen forest. *C. unicolor* was more abundant in the forested plots and *C. taedatus* was more abundant in the tundra plots. This combination allows us to test our model's accuracy for high and low abundance levels in both habitats.


## Species traits

This information was sparse for our selected species - therefore some of the "specific" phenology information is for the family or genus level. 

In general, carabids are ground-dwelling beetles that prefer moist conditions; note that both of our selected species are xerophilous. Most species are omnivorous as larvae and adults, commonly eating both other arthropods and seeds. The life cycle is 1 year long for most species, with adults living 2-4 years. Most taxa have three instars. Pupation is in the ground.  

### *Carabus taedatus* 
<img src="images/CARTAE.png" width="100" />

**Range**: Nearctic; **Habitat**: Boreal forests. Prefers dry, well-drained environments with thin, low vegetation; **Diet**: both larvae and adults are predatory; **Phenology**: probably overwinters as an adult; **Reproduction**: (genus) Females lay their eggs in specially-constructed cells made of mud, twigs and leaves. 

### *Cymindis unicolor*
<img src="images/CYMUNI.png" width="100" />

**Range**: Holarctic; **Habitat**: An arctic-alpine species occurring in tundra and tundra-forest transition zones; prefers dry, treeless, well-drained environments;  **Diet**: no specific information found; **Phenology**: Little known: In Canada, most specimens were found in June and July, but also found as late as October; **Reproduction**: (genus) Females use ovipositors to create cavities in the earth then deposit their eggs. **Remarks**: The species is relatively rare. The *Lebiini* tribe is found in all the major zoogeographical regions of the world; this genus ranges from Costa Rica to arctic tundra.





# Q2: EDA

## Carabid data

### Sampling

Carabid beetles were sampled at 10 plots, with 4 traps per plot. Beetles were collected every two weeks during the summer season (with temperature cutoffs to start and stop sampling). Data were available for 2015--2018 (2019 beetles were not yet identified). 


<img src="output/canopy_height_traps.png" width="2144" />

In 2018, plot 4 was removed from sampling and plot 13 was added. Plot 4 was forested, and plot 13 is a tundra site. Plot 4 was likely removed from sampling due to low number of beetles caught (speculation from Matt Bitters). We took advantage of the new site in 2018 by using plot 13 as a test plot for our model predictions, which Anna will get to later.

### Identification
All individuals were identified to species (or morphospecies) by a parataxonomist. A subset (1511/1974) of those were then sent for expert identification. There was 97% agreement overall between the parataxonomist and the expert taxonomist.
<img src="output/paratax_vs_expertax.png" width="600" />

<!-- ![species_abund_by_year_by_trap](output/paratax_vs_expertax.png) -->


### Species selection 
Most of the species collected were only found rarely, so we first narrowed down the set of species to the seven most abundant (and all were accurately ID'd). 

<img src="output/species_through_years.png" width="3283" />

However, comparing multiple models across seven species is still a bit too much. To narrow further, we selected the two species of these seven that were found both above and below treeline. 

**Selected species**: 

*Carabus taedatus* is more abundant in tundra plots.

<img src="images/CARTAE.png" width="100" />

*Cymindis unicolor* is more abundant in forested plots.

<img src="images/CYMUNI.png" width="50" />

<!-- info/natural history about the species we chose
Do we have an idea for why certain carabids are greater in different nlcdClass given their biology?
-->

We chose these species because we thought it would be interesting to examine model performance for the species that were found in both habitats, especially since they show opposite trends in abundance in the two habitats.

### Abundance over time

<!-- by year -->
<!-- change to just the two species -->
<img src="output/species_abund_by_year_by_trap_nlcd.png" width="2100" />

This plot shows the abundances at each trap in each year for *Carabus taedatus* and *Cymindis unicolor*. *Carabus taedatus* was reliably caught at the tundra sites and was caught on occasion in some forested sites. *Cymindis unicolor* was reliably caught at some (but not all) forested sites and was caught on occasion at one tundra plot (8) and once at another tundra plot (6). 



## Predictor variables  
At the start of this project we spent some time brainstorming and searching the literature for potentially important predictor variables for carabid abundance. Once we had identified a long list of potential predictors, we searched through the available datasets for usable predictors. We ended up using predictor variables from three sources: NEON data included in the carabid data product, NEON AOP data products, and Niwot Ridge LTER environmental data. The predictors we used from each of these sources will be shown below. There were also several predictor variables that we were not able to use: soil moisture (NEON product, not colocated with beetle traps), woody debris (NEON product, not colocated with beetle traps), vegetation (NEON product, difficult to summarize at plot or trap level), and snowmelt date (didn't find or have time to incorporate).

#### Included with carabid data product
Day of year - extracted from collect date  
NLCD class - evergreen forest or grassland herbaceous  
Elevation  

#### AOP (remote sensing) data products 
Data products from NEON's AOP collection were downloaded via the `neonUtilities` package in the form of multiple (over 100) .tif files in a zip folder for each data product. We then merged together all of the .tif's into a single raster in R (`carabids_01_download_AOP.R`) for each data product. These rasters were then used for all of our downstream analyses.

**Canopy height**  
Measured via LIDAR, this data product provides a 1x1m resolution map of canopy height across the site. To summarize for each trap, we took the average canopy height for a 10m radius from each trap location (see`raster::extract` in `carabids_03_EDA_CHM`). 
<img src="output/canopy_height_traps.png" width="600" />



**Leaf area index**  
This is the ratio of leaf area to bare ground area as sensed from above (similar metric to percent canopy cover). This was processed the same way as canopy height, but with some extra data cleaning. Data were available for 2017, 2018, and 2019, but the data looked different year to year. In the end, we used an average of 2017 and 2018 because these issues were difficult to sort out and those years were similar and appeared the most reliable. After contacting someone at NEON, we learned that they used different quality control measures in different years resulting in different values.

**Slope and aspect**  
Slope and aspect were summarized at the plot level (average slope and aspect at 10m radius from plot center). We included them because they may be related to microclimatic factors such as soil moisture (slope) and energy availability via sun exposure (aspect).

<img src="output/slope_traps.png" width="400" />
<img src="output/aspect_traps.png" width="400" />

#### Niwot Ridge LTER data sources
**Precipitation**  
For the forested plots, C1 precipitation data were used; for the tundra plots Saddle precipitation data were used. [NEON also has a precipitation data product, but there were too many missing data points for our use.] We calculated total precipitation during the collection period (2 weeks) for each collect date. This summary gets at the short-term weather patterns that may influence beetle behavior and detectability in traps. [We wanted to get to a long-term precipitation summary, but ran out of time.]
<img src="output/precip_2week_summ.png" width="2240" />

**Growing degree days**  
Using temperature data from the C1 sensor (Saddle sensor data contained missing dates), we calculated growing degree days. Growing degree days were summarized similar to (Nufio et al. 2010) as a cumulative sum of growing degree days in that year until the date of collection. Degree days were quantified for temperatures between 12--33C based on grasshopper physiology at Niwot Ridge (Nufio et al. 2010). We did not have any information regarding the carabid temperature limits, so grasshopper values from the same location were our best estimate.
<img src="output/temp_time.png" width="1600" />

<img src="output/degree_days_cumulative.png" width="1600" />


### Summary of predictors used in models

* Day of year
* NLCD class
* Elevation
* Canopy height (10m radius)
* Leaf area index (LAI) (10m radius)
* Slope
* Aspect
* Precipitation (during trapping period)
* Growing degree days (cumulative over season)



# Q3 and 4: Model building and selection



To build and select a model for abundance of *Carabus taedatus* and *Cymindis unicolor*, our goals were to:  
    1. Identify what are the key predictor variables for each species.   
    2. Find what model best fits the training data species abundance.  

#### Predictor variable selection

We ruled out some variables by evaluating collinearity with the `mctest` package.

![](final_pres_done_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

Note the nearly perfect correlation between growing degree-days and day of year. Comparing all possible variables, we detect overall collinearity between elevation, canopy height, and NLCD class.  


```
## 
## Call:
## imcdiag(x = vars_mat, y = all7sp_dat$sp_abund, method = "VIF")
## 
## 
##  VIF Multicollinearity Diagnostics
## 
##                   VIF detection
## elev          31.2824         1
## LAI_1718avg    2.8618         0
## plot_CHM      17.6503         1
## precip_2weeks  1.0966         0
## DOY            7.4956         0
## plot17aspect   2.4078         0
## plot17slope    1.2610         0
## GDD_cum        6.8395         0
## col_year       1.3224         0
## nlcdclass_num 16.5812         1
## 
## Multicollinearity may be due to elev plot_CHM nlcdclass_num regressors
## 
## 1 --> COLLINEARITY is detected by the test 
## 0 --> COLLINEARITY is not detected by the test
## 
## ===================================
```

After dropping elevation, we failed to detect pairwise collinearity between variables. However, variance inflation factor (VIF) values are above 5 (indicating high correlation) for canopy height, NLCD class, day of year, and cummulative growing degree-day. 


```
## 
## Call:
## imcdiag(x = vars_mat, y = all7sp_dat$sp_abund, method = "VIF")
## 
## 
##  VIF Multicollinearity Diagnostics
## 
##                  VIF detection
## LAI_1718avg   2.1936         0
## plot_CHM      6.3306         0
## precip_2weeks 1.0965         0
## DOY           7.4956         0
## plot17aspect  1.5325         0
## plot17slope   1.0645         0
## GDD_cum       6.8394         0
## col_year      1.3209         0
## nlcdclass_num 9.3633         0
## 
## NOTE:  VIF Method Failed to detect multicollinearity
## 
## 
## 0 --> COLLINEARITY is not detected by the test
## 
## ===================================
```

Canopy height and NLCD class are highly correlated because canopy height is 0 in the tundra. NLCD class is useful in model building since we can vary the intercept by this categorical variable. 

Day of year and growing degree-day are highly correlated. Both capture seasonal variation. 

We don't have enough information right now to select which variable to rule out of each pair. Model fitting with these variables can help us decide. Let's take a step in this direction by utilizing Melissa's `stepGAM` function to compare model AIC values for every iteration of these variables combined. 

Second, use Melissa's stepGAM(), a wrapper function for stepAIC that is compatible with GAMs, to see which predictor variables are important in model fit. 

The effects of plot, trap within plot, canopy height, and slope/aspect are strong for both species. For *Cymindis unicolor* we also see the interaction of day of year and precipitation as an important effect. 

![](final_pres_done_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

#### Model selection

Now that we've identified what are the possible key predictor variables to include in modeling, let's compare a few modeling approaches. We found that GAMMs performed better than GLMMs as they capture well the seasonal flux in carabid abundance. We used the `mgcv` package for fitting GAMMs. Our initial efforts used a poisson or negative binomial distribution to model carabid abundance. However, predictions were severely biased toward 0, so we should next try a zero-inflated poisson model. While we have a number of environmental variables at our disposal to use as predictors, we learned from Harris et al. 2018 that simply an average model may yield the best predictions, even if it's not the best fit model. For this reason, we will compare the model results from the best baseline (variance partitioning) and negative binomial models. First, we iterated through models with many combinations of environmental predictors. 

The best fit model for each species utilized distinct environmental predictor variables. Since we observe clear differences in abundance between NLCD class types for each species, we had the models fit different smooths of continuous variables to each type of NLCD class.  Models for both species used grouping variables that are artifacts of the sampling structure:  
* collection date (24 total)  
* collection year (4 total)  
* plot (10 across site)  
* traps within plots (4 per plot)  

For *Carabus taedatus*, the following environmental predictors were used:  
* day of year  
* canopy height  
* accumulated precipitation over collection period  
* interaction of precipitation and day of year  
* interaction of slope and aspect  

For *Cymindis unicolor*, the following environmental predictors were used:  
* cummulative growing degree days  
* LAI  
* accumulated precipitation over collection period  
* interaction of slope and aspect  
* interaction of precipitation and day of year  
* interaction of growing degree days and precipitation  


```
## sp_abund ~ nlcdClass + s(DOY, bs = "cc", k = 3) + s(precip_2weeks, 
##     by = nlcdClass, bs = "ts", k = 3) + s(plot_CHM, bs = "ts", 
##     k = 3) + te(DOY, precip_2weeks, by = nlcdClass, bs = c("cc", 
##     "ts"), k = 3) + te(plot17aspect, plot17slope, bs = c("ts", 
##     "ts"), k = 3) + s(collectDate, bs = "re") + s(col_year_fac, 
##     bs = "re") + s(plot_trap, bs = "re") + s(plotID, bs = "re")
```

```
## sp_abund ~ nlcdClass + s(GDD_cum, bs = "ts", k = 1) + s(LAI_1718avg, 
##     bs = "ts", k = 1) + s(precip_2weeks, bs = "ts", k = 1) + 
##     te(GDD_cum, precip_2weeks, bs = c("ts", "ts"), k = 3) + te(plot17aspect, 
##     plot17slope, by = nlcdClass, bs = c("ts", "ts"), k = 5) + 
##     te(DOY, precip_2weeks, bs = c("cc", "ts"), k = 4) + s(collectDate, 
##     bs = "re") + s(col_year_fac, bs = "re") + s(plot_trap, bs = "re") + 
##     s(plotID, bs = "re")
```



# Results & Interpretation

## Biological interpretation
Many of the models we selected between had comparable AIC values, so we chose ones that would be most easily interpreted biologically. Now let's interpret them. 

First, let's look at the output of the variance partition model for each species to see where the variability in the data sit.Note that the variance partition model explains over half of the deviance in the data already!


```r
load("data_derived/vp_mod_cartae.Rdata")
load("data_derived/vp_mod_cymuni.Rdata")
summary(vp_mod_cymuni)
```

```
## 
## Family: Negative Binomial(0.448) 
## Link function: log 
## 
## Formula:
## sp_abund ~ s(collectDate, bs = "re") + s(col_year_fac, bs = "re") + 
##     s(plot_trap, bs = "re") + s(plotID, bs = "re")
## 
## Parametric coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  -2.9602     0.6533  -4.531 5.87e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                    edf Ref.df  Chi.sq p-value   
## s(collectDate)  12.014     23  30.763 0.00314 **
## s(col_year_fac)  1.237      3   5.758 0.23388   
## s(plot_trap)    15.193     39  84.973 0.20337   
## s(plotID)        7.473      9 272.018 0.00399 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.164   Deviance explained = 50.9%
## -REML = 359.47  Scale est. = 1         n = 932
```

Now let's look at the best-fit GAMM result plots.


```r
summary(cartae_mod)
```

```
## 
## Family: Negative Binomial(1.107) 
## Link function: log 
## 
## Formula:
## sp_abund ~ nlcdClass + s(DOY, bs = "cc", k = 3) + s(precip_2weeks, 
##     by = nlcdClass, bs = "ts", k = 3) + s(plot_CHM, bs = "ts", 
##     k = 3) + te(DOY, precip_2weeks, by = nlcdClass, bs = c("cc", 
##     "ts"), k = 3) + te(plot17aspect, plot17slope, bs = c("ts", 
##     "ts"), k = 3) + s(collectDate, bs = "re") + s(col_year_fac, 
##     bs = "re") + s(plot_trap, bs = "re") + s(plotID, bs = "re")
## 
## Parametric coefficients:
##                              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                   -3.2960     0.3425  -9.623  < 2e-16 ***
## nlcdClassgrasslandHerbaceous   3.2106     0.4723   6.797 1.07e-11 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                                                          edf Ref.df Chi.sq
## s(DOY)                                             6.947e-05      2   0.00
## s(precip_2weeks):nlcdClassevergreenForest          4.180e-06      2   0.00
## s(precip_2weeks):nlcdClassgrasslandHerbaceous      3.024e-05      2   0.00
## s(plot_CHM)                                        4.338e-06      2   0.00
## te(DOY,precip_2weeks):nlcdClassevergreenForest     1.423e+00      6   3.28
## te(DOY,precip_2weeks):nlcdClassgrasslandHerbaceous 1.820e-04      6   0.00
## te(plot17aspect,plot17slope)                       4.259e-05      8   0.00
## s(collectDate)                                     2.268e-04     23   0.00
## s(col_year_fac)                                    2.356e+00      3  11.18
## s(plot_trap)                                       2.278e-04     38   0.00
## s(plotID)                                          5.161e+00      8  28.57
##                                                     p-value    
## s(DOY)                                              0.40762    
## s(precip_2weeks):nlcdClassevergreenForest           0.40929    
## s(precip_2weeks):nlcdClassgrasslandHerbaceous       0.58624    
## s(plot_CHM)                                         0.60593    
## te(DOY,precip_2weeks):nlcdClassevergreenForest      0.06510 .  
## te(DOY,precip_2weeks):nlcdClassgrasslandHerbaceous  0.42377    
## te(plot17aspect,plot17slope)                        0.49039    
## s(collectDate)                                      0.87587    
## s(col_year_fac)                                     0.00262 ** 
## s(plot_trap)                                        0.98142    
## s(plotID)                                          2.19e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =  0.276   Deviance explained = 51.1%
## -REML = 510.14  Scale est. = 1         n = 932
```

```r
plot(cartae_mod, pages = 5, rug = TRUE, residuals = TRUE, shade=TRUE, shade.col = "lightblue", all.terms = TRUE, seWithMean = TRUE, shift = coef(cartae_mod)[1], scheme=2)
```

![](final_pres_done_files/figure-html/unnamed-chunk-26-1.png)<!-- -->![](final_pres_done_files/figure-html/unnamed-chunk-26-2.png)<!-- -->![](final_pres_done_files/figure-html/unnamed-chunk-26-3.png)<!-- -->


```r
summary(cymuni_mod)
```

```
## 
## Family: Negative Binomial(0.436) 
## Link function: log 
## 
## Formula:
## sp_abund ~ nlcdClass + s(GDD_cum, bs = "ts", k = 1) + s(LAI_1718avg, 
##     bs = "ts", k = 1) + s(precip_2weeks, bs = "ts", k = 1) + 
##     te(GDD_cum, precip_2weeks, bs = c("ts", "ts"), k = 3) + te(plot17aspect, 
##     plot17slope, by = nlcdClass, bs = c("ts", "ts"), k = 5) + 
##     te(DOY, precip_2weeks, bs = c("cc", "ts"), k = 4) + s(collectDate, 
##     bs = "re") + s(col_year_fac, bs = "re") + s(plot_trap, bs = "re") + 
##     s(plotID, bs = "re")
## 
## Parametric coefficients:
##                              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                   -2.6097     0.3951  -6.606 3.95e-11 ***
## nlcdClassgrasslandHerbaceous  -2.1537     0.9249  -2.329   0.0199 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Approximate significance of smooth terms:
##                                                                 edf Ref.df
## s(GDD_cum)                                                7.362e-01      2
## s(LAI_1718avg)                                            2.246e-04      2
## s(precip_2weeks)                                          1.757e-04      2
## te(GDD_cum,precip_2weeks)                                 6.155e-01      4
## te(plot17aspect,plot17slope):nlcdClassevergreenForest     2.899e+00      6
## te(plot17aspect,plot17slope):nlcdClassgrasslandHerbaceous 1.379e+00      2
## te(DOY,precip_2weeks)                                     3.330e+00      9
## s(collectDate)                                            3.628e+00     23
## s(col_year_fac)                                           2.046e+00      3
## s(plot_trap)                                              1.470e+01     38
## s(plotID)                                                 6.563e-05      8
##                                                            Chi.sq  p-value    
## s(GDD_cum)                                                  5.060 0.045303 *  
## s(LAI_1718avg)                                              0.000 0.896302    
## s(precip_2weeks)                                            0.000 0.362979    
## te(GDD_cum,precip_2weeks)                                   1.127 0.126433    
## te(plot17aspect,plot17slope):nlcdClassevergreenForest     193.420 0.002080 ** 
## te(plot17aspect,plot17slope):nlcdClassgrasslandHerbaceous  10.067 0.011047 *  
## te(DOY,precip_2weeks)                                      22.854 0.001168 ** 
## s(collectDate)                                              4.440 0.195747    
## s(col_year_fac)                                             8.322 0.030188 *  
## s(plot_trap)                                               33.590 0.000234 ***
## s(plotID)                                                   0.000 0.946847    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## R-sq.(adj) =   0.15   Deviance explained = 50.1%
## -REML = 351.96  Scale est. = 1         n = 932
```

```r
plot(cymuni_mod, pages = 5, rug = TRUE, residuals = TRUE, shade=TRUE, shade.col = "lightblue", all.terms = TRUE, seWithMean = TRUE, shift = coef(cymuni_mod)[1], scheme=2)
```

![](final_pres_done_files/figure-html/unnamed-chunk-27-1.png)<!-- -->![](final_pres_done_files/figure-html/unnamed-chunk-27-2.png)<!-- -->![](final_pres_done_files/figure-html/unnamed-chunk-27-3.png)<!-- -->


## Prediction
Given the limited size of the carabid dataset over 4 years and at one site, we will validate our models on data from a single plot in a single year. We trained our models on the rest of the data. Remember that in our last sampling year, 2018, NEON swapped out plot 4 for plot 13. We thought it would be interesting to predict beetle abundance at this new plot, given what we know about the other plots and their abundances. 


```r
val_dat <- all7sp_dat %>% 
  filter(plotID == "NIWO_013" & col_year == 2018) 

vp_cartae_pred <- cbind(val_dat %>% filter(para_sciname == "Carabus taedatus") %>% select(sp_abund), data.frame(predict(vp_mod_cartae, newdata = val_dat %>% filter(para_sciname == "Carabus taedatus"), type = "response", se.fit= TRUE)))

vp_cymuni_pred <- cbind(val_dat %>% filter(para_sciname == "Cymindis unicolor") %>% select(sp_abund), data.frame(predict(vp_mod_cymuni, newdata = val_dat %>% filter(para_sciname == "Cymindis unicolor"), type = "response", se.fit= TRUE)))

cartae_pred <- cbind(val_dat %>% filter(para_sciname == "Carabus taedatus") %>% select(sp_abund), data.frame(predict(cartae_mod, newdata = val_dat %>% filter(para_sciname == "Carabus taedatus"), type = "response", se.fit= TRUE)))

cymuni_pred <- cbind(val_dat %>% filter(para_sciname == "Cymindis unicolor") %>% select(sp_abund), data.frame(predict(cymuni_mod, newdata = val_dat %>% filter(para_sciname == "Cymindis unicolor"), type = "response", se.fit= TRUE)))

plot_vp_cartae_pred <- vp_cartae_pred %>%
  ggplot(aes(x=fit, y=sp_abund)) +
  geom_jitter(col="dodgerblue") +
  geom_abline(slope = 1, intercept = 0, col="dodgerblue") +
  ggtitle("Var. Part. model for Carabus taedatus") +
  xlim(0,5) +  ylim(0,5) +
  ylab("Observed count" ) +
  xlab("Predicted count")
plot_vp_cymuni_pred <- vp_cymuni_pred %>%
  ggplot(aes(x=fit, y=sp_abund)) +
  geom_jitter(col="dodgerblue") +
  geom_abline(slope = 1, intercept = 0, col="dodgerblue") +
  ggtitle("Var. Part. model for Cymindis unicolor") +
  xlim(0,5) +  ylim(0,5) +
  ylab("Observed count" ) +
  xlab("Predicted count")
plot_cartae_pred <- cartae_pred %>%
  ggplot(aes(x=fit, y=sp_abund)) +
  geom_jitter(col="dodgerblue") +
  geom_abline(slope = 1, intercept = 0, col="dodgerblue") +
  ggtitle("Best fit model for Carabus taedatus") +
  xlim(0,5) +  ylim(0,5) +
  ylab("Observed count" ) +
  xlab("Predicted count")
plot_cymuni_pred <- cymuni_pred %>%
  ggplot(aes(x=fit, y=sp_abund)) +
  geom_jitter(col="dodgerblue") +
  geom_abline(slope = 1, intercept = 0, col="dodgerblue") +
  ggtitle("Best fit model for Cymindis unicolor") +
  xlim(0,5) +  ylim(0,5) +
  ylab("Observed count" ) +
  xlab("Predicted count")

grid.arrange(plot_vp_cartae_pred,plot_vp_cymuni_pred,plot_cartae_pred,plot_cymuni_pred,nrow=2)
```

![](final_pres_done_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

Our models can predict absence or low abundance, but not more. This is support for trying a zero-inflated poisson next. 


## Cross-validation
Cross-validation is intended to evaluate model performance (esp re: prediction error) using a limited data sample. It estimates the skill of the model on unseen (test) data. There are various cross-validation methods; the key difference is in how they partition the data into training and test datasets. 

We chose the k-fold cross validation method, using *k* = 10. This method has relatively low bias due to random grouping of test and training data. 

### About the method
* This method randomly splits the data into *k* groups for cross-validation

* Then it cycles through the groups, selecting each one in turn as a test data set and using the other *k*-1 groups as a training data set to fit and test the model

* The Mean Square Error [MSE] is calculated and retained for each step

* The mean of all *k* MSEs is reported to evaluate overall model performance 


```r
# all7sp_dat_cartae_k <- all7sp_dat %>% 
#   filter(para_sciname=="Carabus taedatus") %>%
#   mutate(k = sample(1:20, 960, replace=T))
# all7sp_dat_cymuni_k <- all7sp_dat %>% 
#   filter(para_sciname=="Cymindis unicolor" ) %>%
#   mutate(k = sample(1:20, 960, replace=T))
# 
# err_cartae <- rep(NA,960)
# for (i in 1:20 ) {
#     gamcartae_loo <- gam(formula = sp_abund ~ nlcdClass + s(DOY, bs = "cc", k = 3) + 
#                            s(precip_2weeks, by = nlcdClass, bs = "ts", k = 3) + 
#                            s(plot_CHM,bs = "ts", k = 3) + 
#                            te(DOY, precip_2weeks, by = nlcdClass, bs = c("cc", "ts"), k = 3) +
#                            te(plot17aspect, plot17slope, bs = c("ts", "ts"), k = 3) + 
#                            s(collectDate, bs = "re") + s(col_year_fac, bs = "re") + 
#                            s(plot_trap, bs = "re") + s(plotID, bs = "re"), family = nb, 
#                          data = all7sp_dat_cartae_k %>% filter(k != i), method = "REML")
#       
#     err_temp <- ( predict(gamcartae_loo,newdata=all7sp_dat_cartae_k %>% filter(k == i)) -
#                          all7sp_dat_cartae_k %>% filter(k == i) %>% select(sp_abund) )^2 #MSE
#     err_cartae <- rbind(err_cartae, err_temp)
#     rm(gamcartae_loo, err_temp)
# }
# err_cartae <- err_cartae[-1,]
# CV_cartae <- mean(err_cartae)
# 
# err_cymuni <-rep(NA,960)
# for (i in 1:20 ) {
#     gamcymuni_loo <- gam(formula = sp_abund ~ nlcdClass + s(GDD_cum, bs = "ts", k = 1) + 
#                            s(LAI_1718avg, bs = "ts", k = 1) + 
#                            s(precip_2weeks, bs = "ts", k = 1) + 
#                            te(GDD_cum, precip_2weeks, bs = c("ts", "ts"), k = 3) + 
#                            te(plot17aspect, plot17slope, by = nlcdClass, bs = c("ts", "ts"), k = 5) +
#                            te(DOY, precip_2weeks, bs = c("cc", "ts"), k = 4) + 
#                            s(collectDate, bs = "re") + s(col_year_fac, bs = "re") + 
#                            s(plot_trap, bs = "re") + s(plotID, bs = "re"), family = nb, 
#                          data = all7sp_dat_cymuni_k %>% filter(k != i), method = "REML")
#       
#     err_temp <- ( predict(gamcymuni_loo,newdata=all7sp_dat_cymuni_k %>% filter(k == i)) -
#                          all7sp_dat_cymuni_k %>% filter(k == i) %>% select(sp_abund) )^2 #MSE
#     err_cymuni <- rbind(err_cymuni, err_temp)
#     rm(gamcymuni_loo, err_temp)
# }
# err_cymuni <- err_cymuni[-1,]
# CV_cymuni <- mean(err_cymuni)
# 
# err_cartae <- rep(NA,960)
# for (i in 1:20 ) {
#     gamcartae_vp_loo <- gam(formula = sp_abund ~
#                            s(collectDate, bs = "re") + s(col_year_fac, bs = "re") + 
#                            s(plot_trap, bs = "re") + s(plotID, bs = "re"), family = nb, 
#                          data = all7sp_dat_cartae_k %>% filter(k != i), method = "REML")
#       
#     err_temp <- ( predict(gamcartae_vp_loo,newdata=all7sp_dat_cartae_k %>% filter(k == i)) -
#                          all7sp_dat_cartae_k %>% filter(k == i) %>% select(sp_abund) )^2 #MSE
#     err_cartae <- rbind(err_cartae, err_temp)
#     rm(gamcartae_vp_loo, err_temp)
# }
# err_cartae <- err_cartae[-1,]
# CV_cartae_vp <- mean(err_cartae)
# 
# err_cymuni <-rep(NA,960)
# for (i in 1:20 ) {
#     gamcymuni_vp_loo <- gam(formula = sp_abund ~
#                            s(collectDate, bs = "re") + s(col_year_fac, bs = "re") + 
#                            s(plot_trap, bs = "re") + s(plotID, bs = "re"), family = nb, 
#                          data = all7sp_dat_cymuni_k %>% filter(k != i), method = "REML")
#       
#     err_temp <- ( predict(gamcymuni_vp_loo,newdata=all7sp_dat_cymuni_k %>% filter(k == i)) -
#                          all7sp_dat_cymuni_k %>% filter(k == i) %>% select(sp_abund) )^2 #MSE
#     err_cymuni <- rbind(err_cymuni, err_temp)
#     rm(gamcymuni_vp_loo, err_temp)
# }
# err_cymuni <- err_cymuni[-1,]
# CV_cymuni_vp <- mean(err_cymuni)
# 
# # Compare models
# save(CV_cartae, CV_cymuni, CV_cartae_vp, CV_cymuni_vp, file="data_derived/cv_output.Rdata")

load("data_derived/cv_output.Rdata")
CV_mat <- matrix(c(CV_cartae, CV_cymuni, CV_cartae_vp, CV_cymuni_vp), nrow=2, byrow=T)
rownames(CV_mat) <- c("best fit", "var par")
colnames(CV_mat) <- c("Car tae", "Cym uni")
CV_mat
```

```
##            Car tae  Cym uni
## best fit 10.684807 15.25025
## var par   8.724124 12.81825
```


# Discussion - Did we do what we set out to do?

We did meet our primary semester goal:

To create a working model to accurately predict the abundance of 1-2 carabid spp. at NEON's Niwot carabid traps (spatial prediction, Plot 13 in 2018)

This was scaled back quite a bit from our early aspirations. And we did not meet some of our larger goals, including running a working cross-validation

## Future work

If we had more time (and data), we'd include the following items: 

* Once the 2019 carabid data is identified, we'd try to predict our two species for 2019 at all plots and collection dates

* If successful, we'd try to predict all 7 abundant species

* Regardless of outcome, we'd try to use our model to answer the question: "How much data is enough?" for accurate forcasting at different timesteps
