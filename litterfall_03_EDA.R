library(neonUtilities)
library(tidyverse)
library(lubridate)

# clear space
# rm(list = ls())

# load data

# litter <- loadByProduct(dpID="DP1.00003.001", 
#                         site=c("MOAB","ONAQ"),
#                         startdate="2018-05", 
#                         enddate="2018-08")

load(file = "data_derived/litter_woodfall_NIWO.Rdata")

# explore data frames

fielddata <- litter_woodfall$ltr_fielddata
massdata <- litter_woodfall$ltr_massdata
pertrap <- litter_woodfall$ltr_pertrap

head(fielddata)
head(massdata)
head(pertrap)

for ( c in 1:ncol(massdata) ) {
    print(paste("Column name:",names(massdata)[c]))
    print(class(massdata[,c]))
    #print(unique(massdata[,c]))
    #print(range(massdata[,c]))
}

for ( c in 1:ncol(pertrap) ) {
    print(paste("Column name:",names(pertrap)[c]))
    print(class(pertrap[,c]))
    #print(unique(massdata[,c]))
    #print(range(massdata[,c]))
}
unique(massdata$plotID)  
# "NIWO_057" "NIWO_064" "NIWO_058" "NIWO_061" "NIWO_047" "NIWO_062" "NIWO_067" "NIWO_041" 
# "NIWO_051" "NIWO_046" "NIWO_063" "NIWO_040"

range(massdata$setDate)  # "2015-09-26 GMT" "2019-09-25 GMT"
range(massdata$collectDate)  # "2016-06-16 GMT" "2019-10-23 GMT"
unique(year(massdata$collectDate)) # 2016, 2017, 2018, 2019
unique(massdata$collectDate) # ~ monthly, June-November, day/period varied
unique(day(massdata$collectDate))
range(massdata$dryMass) # 0.00 54.23
mean(massdata$dryMass) # 1.180009
median(massdata$dryMass) # 0.01

sum(massdata$dryMass == 0) #537 of 1701
sum(massdata$dryMass >5) #101
sum(between(massdata$dryMass, 0.000000000000001, 1)) #1384 0 and 1, 834 0.00000000000001 and 1

unique(pertrap$plotID) 
#"NIWO_062" "NIWO_051" "NIWO_067" "NIWO_057" "NIWO_041" "NIWO_046" "NIWO_063" "NIWO_064" 
# "NIWO_047" "NIWO_061" "NIWO_058" "NIWO_040"
unique(pertrap$date)  # "2015-07-15 GMT"
unique(pertrap$nlcdClass) # "shrubScrub" "grasslandHerbaceous" "evergreenForest"
unique(pertrap$plotType) #tower
unique(massdata$functionalGroup) 
#  "Leaves", "Woody material", "Seeds", "Needles", "Twigs/branches", "Other", "Flowers", "Mixed"    

# histogram of dry mass collected
ggplot(data=massdata) +
    geom_histogram(mapping = aes(x=dryMass))  

ggplot(data=massdata) +
    geom_boxplot(mapping = aes(y=dryMass, x=plotID))

ggplot(data=massdata) +
    geom_boxplot(mapping = aes(y=dryMass, x=functionalGroup))

## how to get a boxplot (or similar) of the functional group collected at each plot? Not this way...
# ggplot(data=massdata) +
#     geom_boxplot(mapping = aes(y=count(functionalGroup), x=plotID))

# boxplot of dry mass by collection date - wouldn't work, only recognized 1 collection date in plotting although there are several 
ggplot(data=massdata) +
    geom_boxplot(mapping = aes(y=dryMass, x=collectDate))  ## broken....

ggplot(massdata, aes(x=collectDate, y=dryMass, color=plotID)) + 
    geom_point(size=2) + 
    ggtitle("Mass collected by date and plotID")

ggplot(massdata, aes(x=collectDate, y=dryMass, color=functionalGroup)) + 
    geom_point(size=2) + 
    ggtitle("Mass collected by date and veg. group")

ggplot(massdata, aes(x=plotID, y=dryMass, color=functionalGroup)) + 
    geom_point(size=2) + 
    ggtitle("Mass collected by plotID and veg. group")


# density plots weren't useful; and didn't seem to recognize sorting data by collection date or year
