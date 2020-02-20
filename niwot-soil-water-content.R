# Download and visualize soil water content and ion content for Niwot
# this includes 30 minute and 1 minute data at 5 locations

library(neonUtilities)
library(tidyverse)

soil_wc <- loadByProduct("DP1.00094.001", 
                         site = "NIWO", 
                         check.size = FALSE)

# visualize sensor positions (five locations across Niwot)
soil_wc$sensor_positions_00094 %>%
  ggplot(aes(referenceLongitude, referenceLatitude)) + 
  geom_point()

# visualize 30 minute volumetric soil water content (VSWC)
soil_wc$SWS_30_minute %>%
  ggplot(aes(startDateTime, VSWCMean)) + 
  geom_point(size = .1) + 
  facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))

# same for volumetric soil ion content (VSIC)
soil_wc$SWS_30_minute %>%
  ggplot(aes(startDateTime, VSICMean)) + 
  geom_point(size = .1) + 
  facet_wrap(~paste("h:", horizontalPosition, "v:", verticalPosition))
