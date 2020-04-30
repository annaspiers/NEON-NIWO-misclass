# Package ID: knb-lter-nwt.413.10 Cataloging System:https://pasta.edirepository.org.
# Data set title: Air temperature data for Saddle chart recorder, 1981 - ongoing.
# AND 
# Package ID: knb-lter-nwt.411.12 Cataloging System:https://pasta.edirepository.org.
# Data set title: Air temperature data for C1 chart recorder, 1952 - ongoing.
# Data set creator:    - Niwot Ridge LTER 
# Data set creator:  Mark Losleben -  
# Contact:    - Information Manager Niwot Ridge LTER/University of Colorado  - lternwt@colorado.edu
# Contact:  Jennifer Morse -    - jennifer.f.morse@colorado.edu
# Stylesheet v2.7 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 
# Citation: Niwot Ridge LTER and M. Losleben. 2018. Air temperature data for Saddle chart recorder, 1981 - ongoing. ver 10. Environmental Data Initiative. https://doi.org/10.6073/pasta/5cad977568c5d03826f32fdd09c4a069. Accessed 2020-04-25.

# saddle
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-nwt/413/10/afc48b6fab15649c9f91a9367debd2e0" 

# C1
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-nwt/411/12/81101917c12b63474ddef47f104b7128" 

infile1 <- tempfile()
download.file(inUrl1,infile1,method="curl")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                   "LTER_site",     
                   "local_site",     
                   "date",     
                   "airtemp_max",     
                   "flag_airtemp_max",     
                   "airtemp_min",     
                   "flag_airtemp_min",     
                   "airtemp_avg"    ), check.names=TRUE)


# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$LTER_site)!="factor") dt1$LTER_site<- as.factor(dt1$LTER_site)
if (class(dt1$local_site)!="factor") dt1$local_site<- as.factor(dt1$local_site)                                   
# attempting to convert dt1$date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1date<-as.Date(dt1$date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1date) == length(tmp1date[!is.na(tmp1date)])){dt1$date <- tmp1date } else {print("Date conversion failed for dt1$date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1date) 
if (class(dt1$airtemp_max)=="factor") dt1$airtemp_max <-as.numeric(levels(dt1$airtemp_max))[as.integer(dt1$airtemp_max) ]
if (class(dt1$flag_airtemp_max)!="factor") dt1$flag_airtemp_max<- as.factor(dt1$flag_airtemp_max)
if (class(dt1$airtemp_min)=="factor") dt1$airtemp_min <-as.numeric(levels(dt1$airtemp_min))[as.integer(dt1$airtemp_min) ]
if (class(dt1$flag_airtemp_min)!="factor") dt1$flag_airtemp_min<- as.factor(dt1$flag_airtemp_min)
if (class(dt1$airtemp_avg)=="factor") dt1$airtemp_avg <-as.numeric(levels(dt1$airtemp_avg))[as.integer(dt1$airtemp_avg) ]

# Convert Missing Values to NA for non-dates

dt1$airtemp_max <- ifelse((trimws(as.character(dt1$airtemp_max))==trimws("NaN")),NA,dt1$airtemp_max)
dt1$flag_airtemp_max <- as.factor(ifelse((trimws(as.character(dt1$flag_airtemp_max))==trimws("NaN")),NA,as.character(dt1$flag_airtemp_max)))
dt1$airtemp_min <- ifelse((trimws(as.character(dt1$airtemp_min))==trimws("NaN")),NA,dt1$airtemp_min)
dt1$flag_airtemp_min <- as.factor(ifelse((trimws(as.character(dt1$flag_airtemp_min))==trimws("NaN")),NA,as.character(dt1$flag_airtemp_min)))
dt1$airtemp_avg <- ifelse((trimws(as.character(dt1$airtemp_avg))==trimws("NaN")),NA,dt1$airtemp_avg)


# Here is the structure of the input data frame:
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

summary(LTER_site)
summary(local_site)
summary(date)
summary(airtemp_max)
summary(flag_airtemp_max)
summary(airtemp_min)
summary(flag_airtemp_min)
summary(airtemp_avg) 
detach(dt1)               

# write the datafiles to csv for re-use
write.csv(dt1, file="data_derived/Niwot_airtemp_saddle.csv")
write.csv(dt1, file = "data_derived/Niwot_airtemp_C1.csv")
