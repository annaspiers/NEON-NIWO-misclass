

query_api_for_coords = function(strings_to_query, url_string = "http://data.neonscience.org/api/v0/locations/"){
    # takes a vector of strings, appends them to the url_string base, 
    # uses GET to retrieve the data, and organizes the data for each individual string
    
    
    require(dplyr)
    require(httr)
    require(jsonlite)
    
    coords_df = data.frame(locationName = NA, trap.Easting = NA, trap.Northing = NA, coord_uncertainty = NA)
    
    for(i in 1:length(trap_strings)){
        json_list = fromJSON(content(
            GET(url = paste0(url_string, trap_strings[i])), 
            as = "text")) 
        loc_name = json_list$data$locationName
        easting = json_list$data$locationUtmEasting
        northing = json_list$data$locationUtmNorthing
        coord_unc = json_list$data$locationProperties %>% 
            filter(locationPropertyName == "Value for Coordinate uncertainty") %>% 
            pull(locationPropertyValue) %>% 
            as.numeric()
        new_row = data.frame(locationName = loc_name, 
                             trap.Easting = easting, 
                             trap.Northing = northing, 
                             coord_uncertainty = coord_unc)
        coords_df = rbind(coords_df, new_row)
        
    }
    
    return(coords_df[-1,])
}