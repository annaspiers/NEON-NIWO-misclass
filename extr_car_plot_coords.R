extr_car_plot_coords = function(){
    # extracts coordinates for each plot at NIWO and returns them in UTM projection
    
    require(dplyr)
    require(geoNEON)
    
    load("data_raw/carabids_NIWO.Rdata")
    proj = "+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
    
    carabid_spat <- def.extr.geo.os(data = carabid_abund$bet_fielddata, 
                                    'namedLocation', 
                                    locOnly=T) %>%
        dplyr::select(api.decimalLatitude, api.decimalLongitude, Value.for.Plot.ID) %>%
        mutate(api.decimalLatitude = as.numeric(api.decimalLatitude),
               api.decimalLongitude = as.numeric(api.decimalLongitude))
    
    # project data frame to DEM projection
    car_plot_spat = carabid_spat %>% 
        dplyr::select(api.decimalLongitude, api.decimalLatitude) %>%
        SpatialPoints(proj4string = CRS("+proj=longlat")) %>%
        spTransform(proj) 
    
    plot_coords = data.frame('Easting' = car_plot_spat@coords[,1], 'Northing' = car_plot_spat@coords[,2]) %>%
        cbind('plotID' = carabid_spat$Value.for.Plot.ID) %>%
        rename(plot.Easting = Easting, plot.Northing = Northing)
    
    return(plot_coords)
}