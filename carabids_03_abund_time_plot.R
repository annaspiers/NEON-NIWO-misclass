library(dplyr)
library(ggplot2)

clean_dat <- read.csv("data_derived/model_df_by_species_in_sample.csv")

# Plot taxa's abundance at each trap through time

clean_dat %>% 
    filter(para_sciname %in% c('Cymindis unicolor','Carabus taedatus')) %>%
    group_by(para_sciname, plot_trap, col_year, sp_abund, nlcdClass) %>%
    summarize(sum_abund = sum(sp_abund)) %>%
    mutate(occ = as.character(ifelse(sum_abund>0,1,0))) %>%
    ggplot(aes(x = col_year, 
               y = plot_trap, 
               color = nlcdClass)) +
    geom_point(aes(size = sum_abund)) +
    scale_size_area() + #value 0 will have zero size
    scale_color_manual(values=c("darkgreen","lightblue")) +
    theme_bw() +
    theme(axis.text.x = element_text(size=6)) +
    facet_grid(. ~ para_sciname)
ggsave("output/species_abund_by_year_by_trap_nlcd.png", width = 7, height = 6)