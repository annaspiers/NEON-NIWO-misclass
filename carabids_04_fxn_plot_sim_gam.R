plot_sim_gam = function(sims, mod_obj, obs_col = 'sp_abund'){
    require(ggplot2)
    require(dplyr)
    
    fitted_obs = data.frame('observed' = mod_obj$data[[obs_col]], 
                            'fitted_vals' = mod_obj$fitted.values)
    
    p = as.data.frame(sims) %>% 
        mutate(observed = mod_obj$data[[obs_col]]) %>%
        melt(id = 'observed') %>% 
        rename(simulated = 'value') %>% 
        ggplot() +
        geom_jitter(aes(x = simulated, y = observed), color = 'gray') +
        geom_point(data = fitted_obs, aes(x = fitted_vals, y = observed)) +
        xlab('fitted values (black), simulated (gray)') +
        theme_bw()
    return(p)
}