# This script takes the cleaned NEON carabid data as input and returns variables
# for JAGS data and initialization. This script includes output for the full and
# reduced models 

library(tidyr) #uncount()
library(tibble) #rownames_to_column()
library(dplyr)

return_jags_input <- function(model, alpha_nondiag = 2, alpha_diag = 200) {
    
    # Load NEON Niwot Ridge carabid data 
    all_paratax_by_ind <- readRDS("data/all_paratax_df.rds") %>%
        rename(para_morph_combo = scimorph_combo)  %>%
        uncount(individualCount) %>%
        rownames_to_column() 
    all_paratax_by_ind$exp_sciname <- NA
    expert_df <- readRDS("data/expert_df.rds") 
    pinned_df <- readRDS("data/pinned_df.rds") 
    expert_pinned_df <- expert_df %>%
        rename(exp_sciname = scientificName) %>%
        left_join(pinned_df %>% 
                      dplyr::select(individualID, subsampleID, 
                                    para_morph_combo = scimorph_combo))
    rm(expert_df, pinned_df)
    
    # List 1) the unique expert taxonomist ID's and 2) the parataxonomist IDs that the
    # expert taxonomist didn't use. 
    # These identifications will be the row names of theta and M
    rownames <- c(expert_pinned_df %>%   #1) the unique expert taxonomist ID's
                      distinct(exp_sciname) %>%
                      pull(exp_sciname),
                  all_paratax_by_ind %>%   #2) the parataxonomist IDs that the expert taxonomist didn't use. 
                      distinct(scientificName) %>%
                      filter(!scientificName %in% unique(expert_pinned_df$exp_sciname)) %>% 
                      pull(scientificName) ) 
    rownames <- sort(rownames)
    
    # List 1) the unique parataxonomist and morphospecies ID's and 2) the expert taxonomist IDs that the
    # parataxonomist didn't use. 
    # These identifications will be the column names of theta and M
    extra_colnames <- c(all_paratax_by_ind %>%   #1) the unique parataxonomist and morphospecies ID's
                            distinct(para_morph_combo) %>%
                            pull(para_morph_combo),
                        expert_pinned_df %>%   #2) the expert taxonomist IDs that the parataxonomist didn't use. 
                            distinct(exp_sciname) %>%
                            filter(!exp_sciname %in% unique(all_paratax_by_ind$para_morph_combo)) %>% 
                            pull(exp_sciname) ) 
    extra_colnames <- sort(extra_colnames)
    extra_indices <- which(!extra_colnames %in% rownames)
    colnames <- c(rownames, extra_colnames[extra_indices])
    rm(extra_colnames, extra_indices)
    
    assertthat::assert_that(all(rownames == colnames[1:length(rownames)]))
    
    # Join parataxonomist and expert ID tables ---------------------------
    # Ex: If parataxonomist counts 5 animals in one subsample & expert IDs two
    # individuals, assign known species IDs to 2 of the 5 rows at random (this is
    # valid because the 5 individuals are exchangeable)
    
    # initialize counter for for-loop
    subsamps <- expert_pinned_df %>%
        distinct(subsampleID) %>% pull(subsampleID)
    # para_new will replace the all_paratax_df since dplyr doesn't play nice in for-loops
    # (I couldn't assign to a filtered object)
    # initialize para_new with subsamples that the expert ID never looked at
    para_new <- all_paratax_by_ind %>% 
        filter(!(subsampleID %in% subsamps)) 
    for (id in subsamps) {
        temp_exp_df <- expert_pinned_df %>%
            dplyr::select(subsampleID, exp_sciname) %>%
            filter(!is.na(exp_sciname), subsampleID == id)
        temp_para_df <- all_paratax_by_ind %>% 
            filter(subsampleID == id)
        for (row in 1:nrow(temp_exp_df)) {
            # check whether exp_sciname is empty in first row of paratax-df
            while (!is.na(temp_para_df$exp_sciname[row])) {
                row = row + 1
            } 
            # assign exp_sciname from pinned-expert-combo-df to corresponding row of paratax-df
            temp_para_df$exp_sciname[row] <- temp_exp_df$exp_sciname[row]
        }
        para_new <- rbind(para_new, temp_para_df)
    }
    para_new <- para_new %>% as_tibble()
    
    ## Sanity check
    assertthat::assert_that(sum(!is.na(para_new$exp_sciname)) == nrow(expert_pinned_df))
    # Choose a subsample that has fewer expert IDs than sorting IDs: ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=
    expert_pinned_df %>%
        filter(!is.na(exp_sciname), subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
    para_new %>% filter(subsampleID == "ta53RQUlpj5WhQksfJeD+QBAxMj6BQGBllkC8fUqt68=")
    # We should see 6 out of 7 individuals in para_new have an assigned exp_sciname
    
    # Define alpha
    alpha <- matrix(alpha_nondiag, nrow = length(rownames), ncol = length(colnames))
    diag(alpha) <- alpha_diag
    
    if (model == "reduced") {
        # Keep only expert IDed individuals for reduced model
        para_new <- para_new %>%
            filter(!is.na(exp_sciname))
    }
    
    # Define L
    L <- reshape2::acast(para_new, plotID ~ collectDate ~ col_year)
    
    if (model == "full"){
        
        ## "Ground truth" data
        # We have a subset of the data with known species IDs from expert identification
        # We partly observe k, the expertID column in y_df
        # We partly observe z
        z.dat <- array(NA, dim = c(length(unique(para_new$plotID)),
                                   length(rownames),
                                   length(unique(all_paratax_by_ind$col_year))),
                       dimnames = list(sort(unique(para_new$plotID)), #plots
                                       rownames, #expert IDs
                                       unique(all_paratax_by_ind$col_year)))
        # Grab values from casted z.dat array and fill in values in final z.dat array
        z.dat_cast <- expert_pinned_df %>%
            mutate(occ = 1) %>%
            reshape2::acast(plotID ~ exp_sciname ~ col_year,
                            fill=-999, drop=F, value.var = "occ")
        z.dat_cast[z.dat_cast == -999] <- NA
        z.dat_cast[z.dat_cast > 0] <- 1
        
        for(plot in dimnames(z.dat_cast)[[1]]) {
            for(spec in dimnames(z.dat_cast)[[2]]) {
                for(year in dimnames(z.dat_cast)[[3]]) {
                    z.dat[plot,spec,year] <- z.dat_cast[plot,spec,year]
                }
            }
        }
        rm(z.dat_cast,plot,spec)
        
        # Initialize Z
        z.init <- z.dat
        for (i in 1:dim(z.init)[1]) {
            for (t in 1:dim(z.init)[3]) {
                z.init[i,,t] <- sample(c(0,1), replace=TRUE, size=dim(z.init)[2])
            }
        }
        # initialize known values as NA, otherwise model will throw error
        z.init[z.dat == 1] <- NA
        
        # Check that where L>0 for a species, z.init>0 for that species/site/year combo
        for (i in 1:dim(z.init)[1]) {
            for (k in 1:dim(z.init)[2]) {
                for (t in 1:dim(z.init)[3]) {
                    if (sum(L[i,,t], na.rm = TRUE) > 0 ) {
                        ifelse(z.init[i,k,t] == 0, 1, z.init[i,k,t])
                    }
                }
            }
        }
        rm(i,k,t)
    }
    
    ## Imperfect species classifications 
    # Probability vector `y` for with a record for each detection. 
    # Noisy classifier with skill that vary by species.
    y_df <- para_new %>%
        dplyr::select(plotID, collectDate, parataxID = para_morph_combo, 
                      expertID = exp_sciname, col_year)
    y_df <- y_df %>%
        left_join(y_df %>%
                      group_by(plotID) %>%
                      summarize(n=n()) %>%
                      mutate(plotID_idx = 1:n()) %>%
                      dplyr::select(plotID, plotID_idx)) %>%
        mutate(parataxID_idx = match(y_df$parataxID, colnames),
               expertID_idx = match(y_df$expertID, rownames)) 
    
    
    if (model=="full") {
        inputlist <- list(rownames=rownames, colnames=colnames, para_new=para_new, 
                          alpha=alpha, L=L, y_df=y_df, z.dat=z.dat, z.init=z.init)
        return(inputlist)
    }
    
    if (model=="reduced") {
        inputlist <- list(rownames=rownames, colnames=colnames,
                          alpha=alpha, L=L, y_df=y_df)
        return(inputlist)
    }
    
}




