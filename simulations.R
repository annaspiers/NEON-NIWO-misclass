# Simulation varying fraction of samples validated 

# This simulation demonstrates the model behavior when varying fractions of
# samples are validated at the individual level. Data are generated using model
# priors. We implement the simulations using a single-season version of the featured
# classification-occupancy model. This simulation illustrates model behavior on
# a generic dataset where all samples (e.g., individuals, states) are
# imperfectly classified so as to reach a broader audience of ecologists that
# are generally interested in classification problems, regardless of
# application. 

# AIS revisit the following to make sure we actually do these:
# We show how this model improves ecological estimation/inference by
# illustrating how close the estimated occupancy or demographic rates are to the
# true parameter values. We also show what fraction of data need to be validated to
# achieve certain model accuracy.

# Model assumptions
# Occupancy has closure throughout the season

library(jagsUI)
library(MCMCvis)
library(ggmcmc)
library(dplyr)
library(tibble) #rownames_to_column
library(parallel)
library(MASS) #mvrnorm()
library(boot) #inv.logit()
library(MCMCpack) #rdirichlet()
library(extraDistr) #rcat()
options(digits = 3)

# library(data.table)
# temp <- list.files(path = "output/simulations", pattern = "sim_out_*", full.name=TRUE)
# sims_jm_list <- lapply(temp, readRDS) #takes 1 min
# sims_jm <- rbindlist(sims_jm_list)
# remaining <- setdiff(1:max(sims_jm$dataset),sims_jm$dataset) 
# rm(temp, sims_jm_list, sims_jm)


# Simulate data from model priors
# By generating data from the model, fitting the model to these data, then computing posteriors for (occupancy, demographic rate) parameter estimates, we have theoretical guarantees about interval coverage. 

# Generate data from model priors --------------------------------------------------

sim_1dataset <- function(iter) {
    nsite <- 30 #number of sites
    nsurv <- 3  #number of surveys in the season
    K_imp <- 3  #number of states identified by imperfect classifier
    K_val <- 2   #number of states identified in validation
  
    ### Generate parameter values from model priors
    mu_spec <- rnorm(2, mean=0, sd=1)
    mu_site <- rep(0,2)
  
    R <- diag(rep(1, 2)) # 2 dim since we're estimating psi and lambda
    Tau_spec <- Tau_site <- rWishart(n=1, df=10, Sigma=R)
    eps_spec <- mvrnorm(n=K_val, mu=mu_spec, Sigma=solve(Tau_spec[,,1], R) )
    eps_site <- mvrnorm(n=nsite, mu=mu_site, Sigma=solve(Tau_site[,,1], R) )
  
    # Psi (occupancy probability), Lambda (encounter rate)
    psi <- lambda <- matrix(NA,nrow=nsite,ncol=K_val)
    for (i in 1:nsite){
      for (k in 1:K_val) {
        psi[i, k] <- inv.logit(eps_site[i, 1] + eps_spec[k, 1])
        lambda[i, k] <- exp(eps_site[i, 2] + eps_spec[k, 2])
      }
    }
  
    # Theta (classification probability)
    alpha <- matrix(100, nrow=K_val, ncol=K_imp)
    diag(alpha) <- 600 # couldn't find great match with 2/3 species here to
    Theta <- matrix(NA, nrow=K_val, ncol=K_imp)
    for (k in 1:K_val) {
      Theta[k,1:K_imp] <- rdirichlet(1,alpha[k,1:K_imp])
    }
  
    ### Simulate datasets for classification-occupancy model
    # Z (expected occupancy) and zlam (expected abundance given occupancy)
    z <- zlam <- matrix(NA, nrow=nsite, ncol=K_val)
    for (i in 1:nsite) {
      for (k in 1:K_val) {
        z[i, k] <- rbinom(1, 1, psi[i, k]) #bernoulli trial
        zlam[i, k] <- z[i, k] * lambda[i, k]
      }
    }
  
    # L (expected number of encountered individuals)
    L <- matrix(NA, nrow=nsite, ncol=nsurv)
    for (i in 1:nsite) {
      for (j in 1:nsurv) {
        L[i, j] <- rpois(1, sum(zlam[i, 1:K_val]) ) #sum zlam across each row for number of imperfectly detected individuals
      }
    }
    Ltot <- sum(L)

    # Pi (dirichlet distribution for each encountered individual)
    # k (validated species for an individual)
    # y (imperfectly observed state for an individual)
    site <- sample(x=which(rowSums(z)>0), size=Ltot, replace=T, prob=rowSums(zlam)[rowSums(zlam)!=0]) #sample(1:nsite, Ltot, replace=TRUE)
    pi <- matrix(NA, nrow=Ltot, ncol=K_val)
    k <- y <- rep(NA, Ltot)
    for (l in 1:Ltot) {
      pi[l, 1:K_val] <- zlam[site[l], 1:K_val] / sum(zlam[site[l], 1:K_val])
      k[l] <- rcat(1, pi[l, 1:K_val])
      y[l] <- rcat(1, Theta[k[l], 1:K_imp])
    }

    # Vary the fraction of individuals that are validated ---------------------

    # For this simulation, we are changing the fraction of samples that were validated.
    # At this point, the object k stores the validated state value for 100% of samples.
    # We need to vary the the number of values in k that are NA, thereby varying the fraction
    # of individuals that are validated
    n <- 100
    fracs <- sort(runif(n, min=0.01, max=0.99), decreasing=T) 
    
    k_mat <- matrix(data=k, nrow=Ltot, ncol=length(fracs))
    colnames(k_mat) <- paste0("k",fracs)
    for (j in 1:ncol(k_mat)) {
      if (j!=1) {k_mat[,j] <- k_mat[,j-1]}
      alreadyNA <- sum(is.na(k_mat[,j]))
      needtobeNA <- round((1-fracs[j])*Ltot)
      ind <- sample(which(!is.na(k_mat[,j])), (needtobeNA - alreadyNA), replace=F)
      k_mat[ind,j]<- NA
    }
    
    val_perc <- 0.1 #fraction of samples to select to calculate validation metrics
    # Prepare input data for calculating validation metrics
    # To validate models, withhold imperfect IDs from individuals that have been validated
    #perc_withld <- 0.2 #percent individuals that have been validated that will have imperfect ID withheld
    num_withld <- ifelse(round(val_perc*Ltot)==0,1,round(val_perc*Ltot)) #take 10% of the total dataset to have imperfect IDs set to NA
    # Select first fracs index that has at least num_withld NAs
    first_frac_for_val <- min(which(round((1-fracs)*Ltot) >= num_withld,arr.ind=TRUE)) #index of first fraction < 1-frac_val
    # Confirm 
    assertthat::assert_that(sum(is.na(k_mat[,first_frac_for_val])) >= num_withld)
    assertthat::assert_that(sum(is.na(k_mat[,first_frac_for_val-1])) < num_withld)
    
    # Select rows at random that will have paratax ID withheld
    rows_withld <- data.frame(k=k_mat[,first_frac_for_val]) %>%
      rownames_to_column() %>%
      filter(is.na(k)) %>%
      slice_sample(n=num_withld) %>% #prop=perc_withld) %>%
      pull(rowname) #select vector of row indices that will be withheld

   # Initialize df to store jmsummm for entire dataset (includes all k values)
    jmsumm_iter <- data.frame()
    val_met_iter <- data.frame()

  for (current in 1:length(fracs)) {
    op <- options(digits.secs = 2)
    start_time <- Sys.time()

    # Specify which k_frac will be used to initialize z.init, z.dat, and run the model
    # manually specify which column index of k_mat
    k_sim <- k_mat[,current]
    
    # Partially observed occupancy --------------------------------------------
    z.dat_df <- data.frame(site,k_sim) %>%
      filter(!is.na(k_sim)) #filter out the NAs so that df contains only validated individuals
    z.dat_df$occ <- 1
    # cast into a nsite x k_val matrix that is full of NAs and 1s
    z.dat_cast <- z.dat_df %>%
      reshape2::acast(site ~ k_sim, fill=-999, drop=F, value.var = "occ")
    z.dat_cast[z.dat_cast == -999] <- NA
    z.dat_cast[z.dat_cast > 0] <- 1

    # Initialize z.dat to have full dimensions nsite x K_val. A column/row can be lost if all values are NA
    z.dat <- array(NA, dim = c(nsite, K_val), dimnames = list(1:nsite, 1:K_val))
    for (plot in dimnames(z.dat_cast)[[1]]) {
      for (spec in dimnames(z.dat_cast)[[2]]) {
        z.dat[plot,spec] <- z.dat_cast[plot,spec]
      }
    }

    # Initialize Z
    z.init <- z
    z.init[z.dat == 1] <- NA # initialize known values as NA, otherwise model will throw error

    # Check that where L>0 for a species, z.init>0 for that species/site/year combo
    for (i in 1:dim(z.init)[1]) {
      for (j in 1:dim(z.init)[2]) {
        if (sum(L[i,], na.rm = TRUE) > 0 ) {
          ifelse(z.init[i,j] == 0, 1, z.init[i,j])
        }
      }
    }

    # Fit model ---------------------------------------------------------------
    nc <- 4
    ni <- 1200 #1000
    na <- 400 #200

    # Run full model
    jags_d_full <- list(nsite = dim(L)[1],
                        nsurv = dim(L)[2],
                        K_val = dim(alpha)[1],
                        K_imp = dim(alpha)[2],
                        L = L,
                        alpha = alpha,
                        Ltot = sum(L),
                        site = site,
                        k = k_mat[ ,current], #true ID is known if  individual was labeled by expert
                        y = y, #for all individuals, we get paratxonomist IDs
                        z = z.dat,
                        R = R)
    jm_full <- jags(data = jags_d_full,
                    inits = function(){ list(z = z.init) },
                    parameters.to.save = c("psi","lambda", "Theta", "k"),
                    model.file = "sim_full_JAGS.txt",
                    n.chains = nc,
                    n.adapt = na,
                    n.iter= ni,
                    n.burnin = na)
    jmsumm_full <- MCMCsummary(jm_full) %>%
      mutate(fraction = fracs[current],
             dataset = iter,
             model = "full") %>%
      rownames_to_column()

    # Run reduced model
    jags_d_red <- list(K_val = dim(alpha)[1],
                       K_imp = dim(alpha)[2],
                       alpha = alpha,
                       k = data.frame(y, k_sim) %>%
                         filter(!is.na(k_sim)) %>% 
                         pull(k_sim), #true ID is known if  individual was labeled by expert
                       y = data.frame(y=y, k=k_sim) %>%
                         filter(!is.na(k_sim)) %>% 
                         pull(y), #for all individuals, we get paratxonomist IDs
                       Ltot = nrow(data.frame(k_sim) %>% filter(!is.na(k_sim))))
    jm_red <- jags(data = jags_d_red,
                   inits = NULL,
                   parameters.to.save = c("Theta"), #no psi or lambda to estimate with full
                   model.file = "sim_red_JAGS.txt",
                   n.chains = nc,
                   n.adapt = na,
                   n.iter= ni,
                   n.burnin = na)
    jmsumm_red <- MCMCsummary(jm_red) %>%
      mutate(fraction = fracs[current],
             dataset = iter,
             model = "reduced") %>%
      rownames_to_column()

    # Initialize dataframes that store true parameter values
    psi_df <- data.frame(param = "psi",
                         i=rep(1:nsite,K_val),
                         j=rep(1:K_val,each=nsite),
                         true = c(psi),
                         dataset = iter)
    lambda_df <- data.frame(param = "lambda",
                            i=rep(1:nsite,K_val),
                            j=rep(1:K_val,each=nsite),
                            true = c(lambda),
                            dataset = iter)
    Theta_df <- data.frame(param = "Theta",
                           i=rep(1:K_val,K_imp),
                           j=rep(1:K_imp,each=K_val),
                           true = c(Theta),
                           dataset = iter)
    trueparams_df <- rbind(psi_df, lambda_df, Theta_df)

    est_true_df_full <- jmsumm_full %>%
      as_tibble %>%
      filter(rowname != "deviance") %>%
      mutate(param = sub("\\[.*","",rowname)) %>%
      relocate(param) %>%
      filter(param!="k") %>%
      separate("rowname", into = c("i", "j"), sep = ",") %>%
      mutate_at(c('i', 'j'), readr::parse_number) %>%
      left_join(trueparams_df)

    est_true_df_red <- jmsumm_red %>%
      as_tibble %>%
      filter(rowname != "deviance") %>%
      mutate(param = sub("\\[.*","",rowname)) %>%
      relocate(param) %>%
      separate("rowname", into = c("i", "j"), sep = ",") %>%
      mutate_at(c('i', 'j'), readr::parse_number) %>%
      left_join(trueparams_df)

    # Append df's together
    jmsumm_iter <- rbind(jmsumm_iter, est_true_df_full, est_true_df_red)
    
    
    ### Validation metrics
    if (current >= first_frac_for_val) {
      # Process validation model output for metrics
      # Filter to validation data
      holdout <- data.frame(y,k,k_sim) %>% 
        rownames_to_column() %>%
        rename(full_individual=rowname) %>%
        mutate(val = ifelse(full_individual %in% rows_withld, 1, 0)) %>%
        filter(val==1) %>% 
        dplyr::select(full_individual,y,real_expID=k) 
      
      # Posterior draws of predicted classifications
      full_k_out <- MCMCchains(jm_full, params='k') %>% #y to k
        as_tibble() %>%
        dplyr::select(matches(paste0("k\\[",holdout$full_individual,"\\]"))) %>% #y to k
        mutate(draw=1:n()) %>%
        pivot_longer(cols=-draw,names_to="full_individual",values_to="est_expID") %>% #changed from "predID" to "est_impID" to est_expID
        mutate_at(c('full_individual'), readr::parse_number) %>%
        mutate(full_individual=as.character(full_individual)) %>%
        left_join(holdout %>% dplyr::select(full_individual,real_expID)) %>% #added real_impID
        mutate(match = (est_expID==real_expID),#match = (predID==trueID),
               model = "full")

      # Generate confusion matrices for each posterior draw ---------------------
      #max_draw <- 2000 #what does this do?
      all_combos <- expand.grid(draw = 1:max(full_k_out$draw),
                                est_impID = 1:max(y), #ais changed predID_idx to predID, same from trueID
                                trueID = 1:max(y))  #%>% as_tibble %>% filter(draw < max_draw)

      full_cm <- full_k_out %>%
        #filter(draw < max_draw) %>%
        count(draw, est_expID, real_expID) %>% #predID to est_impID
        full_join(all_combos) %>% # fill in implicit zeros
        reshape2::acast(draw ~ real_expID ~ est_expID, #predID to est_impID
                        value.var = "n",
                        fill = 0)
      assertthat::assert_that(!any(is.na(full_cm)))                # no NA vals
      assertthat::assert_that(dim(full_cm)[2] == dim(full_cm)[3])  # square matrices

      get_metrics <- function(confusion_matrix) {
        # confusion_matrix is a (true, pred) square matrix
        true_positives <- diag(confusion_matrix)
        false_positives <- colSums(confusion_matrix) - diag(confusion_matrix)
        false_negatives <- rowSums(confusion_matrix) - diag(confusion_matrix)
        precision <- true_positives / (true_positives + false_positives)
        recall <- true_positives / (true_positives + false_negatives)
        f1 <- 2 * (precision * recall) / (precision + recall)
        tibble(ID = 1:length(f1),
               precision = precision,
               recall = recall,
               f1 = f1)
      }

      full_metrics <- apply(full_cm, 1, get_metrics) %>%
        bind_rows(.id = "draw") %>%
        mutate(model = "full")

      # Create code to list accuracy and macro-averages in one line
      val_metrics <- full_metrics %>%
        group_by(model,ID) %>%
        #AIS one thing to look into for why the validaiton metrics look whack is whether I should be taking the average
        #of a single species' metrics across draws before calculating the macro-average for the model
        summarize(precision=mean(precision, na.rm=T), recall=mean(recall, na.rm=T), f1=mean(f1, na.rm=T)) %>%
        #filter(ID != 3) %>%
        #group_by(model) %>%
        summarize(precision=mean(precision, na.rm=T), recall=mean(recall, na.rm=T), f1=mean(f1, na.rm=T)) %>%
        left_join(full_k_out %>% #full_join(full_k_out,reduced_k_out) %>%
                    group_by(model) %>%
                    summarize(accuracy=mean(match))) %>%
        mutate(fraction = fracs[current],
               dataset = iter) %>%
        relocate(accuracy, .after=model) %>%
        relocate(model, .after=dataset)
      
      val_met_iter <- rbind(val_met_iter, val_metrics)
    }
    
    # Write cluster status to a txt file
    write(paste0("Finished dataset ",iter," with ",fracs[current]," samples validated. ",current,"/",n," fractions run for this dataset. ",Sys.time()),file="monitoring.txt",append=TRUE)
  }
  # Save iter-specific df's locally
  saveRDS(jmsumm_iter, paste0("output/simulations/sim_out_",iter,".rds")) #true params and estimates 
  saveRDS(val_met_iter, paste0("output/simulations/sim_val_met_",iter,".rds")) #validation metrics 
}

# Parallelize processing
cl <- makeCluster(12, type = "PSOCK") #12 is better than 24
cl #Print cluster info
clusterEvalQ(cl, c(library(jagsUI), # Send libraries to clusters
                   library(MCMCvis),
                   library(ggmcmc),
                   library(dplyr),
                   library(tibble), 
                   library(parallel),
                   library(MASS), 
                   library(boot), 
                   library(MCMCpack),
                   library(extraDistr)) )
clusterEvalQ(cl, "sim_full_JAGS.txt") # Send jags txt to clusters
clusterEvalQ(cl, "sim_red_JAGS.txt") 

# Run simulations
set.seed(201023923)
iter_tot <- 1000000
system.time( 
  clusterApplyLB(cl, 
               x= c(remaining), #5620:iter_tot,  #setdiff(1:max(sims_jm$dataset),sims_jm$dataset), 
               sim_1dataset)
) #use top to monitor progress

stopCluster(cl)




# SCRAP -------------------------------------------------------------------

# goes under 'assertthat' lines
# This code was needed due to a mistake I made in mis-naming two of the variables in all_combos
# if (dim(full_cm)[2] != dim(full_cm)[3]) {
#   full_cm_temp <- array(NA, dim = c(dim(full_cm)[1],3,3)) #initialize
#   reduced_cm_temp <- array(NA, dim = c(dim(full_cm)[1],3,3)) #initialize
#   if (dim(full_cm)[2]==1) {
#     for (i in 1:dim(full_cm)[3]) {
#       # Add a column for the 3rd imperfect ID
#       full_cm_temp[,,i] <-    cbind(full_cm[,,i], rep(0,dim(full_cm)[1]), rep(0,dim(full_cm)[1]))
#       reduced_cm_temp[,,i] <- cbind(reduced_cm[,,i], rep(0,dim(full_cm)[1]), rep(0,dim(full_cm)[1]))
#     }
#   } else if (dim(full_cm)[2]==2) {
#     for (i in 1:dim(full_cm)[3]) {
#       # Add a column for the 3rd imperfect ID
#       full_cm_temp[,,i] <-    cbind(full_cm[,,i], rep(0,dim(full_cm)[1])) 
#       reduced_cm_temp[,,i] <- cbind(reduced_cm[,,i], rep(0,dim(full_cm)[1]))
#     }
#   }
#   full_cm <- full_cm_temp
#   reduced_cm <- reduced_cm_temp
# }
# assertthat::assert_that(dim(full_cm)[2] == dim(full_cm)[3])  # square matrices


# n <- 100
# #fracs <- runif(n, min=0.01, max=0.99), decreasing=T) 
# fracs <- numeric() 
# # The following for loop is to prevent 0 individuals from being withheld for validation.
# # If frac is too small, then taking 10% of that fraction of the individuals in the dataset will yield 0 rows.
# # Prevent this by generating a new frac if this scenario arises.
# for (i in 1:n) {
#   frac <- runif(1, min=0.1, max=0.99) #min=0.01
#   while (floor(frac*0.1*Ltot)<2) {
#     frac <- frac+.01
#     # Rather than re-sampling a new frac, add 0.01 to the frac until needed so we don't create a bias with not enough low fracs
#   }
#   fracs <- c(fracs, frac)
# }
# fracs <- sort(fracs, decreasing=T)



# ## Run models for validation metrics
# # Full model
# jags_d_full_val <- list(nsite = dim(L)[1],
#                     nsurv = dim(L)[2],
#                     K_val = dim(alpha)[1],
#                     K_imp = dim(alpha)[2],
#                     L = L,
#                     alpha = alpha,
#                     Ltot = sum(L),
#                     site = site,
#                     k = k_df_val$expID_val, #k_mat[ ,current],
#                     y = y, #y_df_val$impID_val, 
#                     z = z.dat,
#                     R = R)
# jm_full_val <- jags(data = jags_d_full_val,
#                 inits = function(){ list(z = z.init)},
#                 parameters.to.save = c("Theta", "k"), #y to k
#                 model.file = "sim_full_JAGS.txt",
#                 n.chains = nc,
#                 n.adapt = na,
#                 n.iter= ni,
#                 n.burnin = na)
# #MCMCtrace(jm_full_val) 
# 
# # Reduced model
# jags_d_red_val <- list(K_val = dim(alpha)[1],
#                    K_imp = dim(alpha)[2],
#                    alpha = alpha,
#                    
#                    # k = k_df_val %>%
#                    #   filter(!is.na(k_sim)) %>% #k_sim to expID_val
#                    #   pull(expID_val), #k_sim to expID_val
#                    # y = k_df_val %>%
#                    #   filter(!is.na(k_sim)) %>% #k_sim to expID_val
#                    #   pull(y), #impID_val to y
#                    # Ltot = nrow(k_df_val %>% filter(!is.na(k_sim)))) #k_sim to expID_val
#                    k = k_df_val %>%
#                      filter(!is.na(expID_val)) %>% #k_sim to expID_val
#                      pull(expID_val), #k_sim to expID_val
#                    y = k_df_val %>%
#                      filter(!is.na(expID_val)) %>% #k_sim to expID_val
#                      pull(y), #impID_val to y
#                    Ltot = nrow(k_df_val %>% filter(!is.na(expID_val)))) #k_sim to expID_val
# jm_red_val <- jags(data = jags_d_red_val, 
#                inits = NULL, #function(){ list(k = KINIT)}, #NULL
#                parameters.to.save = c("Theta", "k"), #y to k
#                model.file = "sim_red_JAGS.txt",
#                n.chains = nc,
#                n.adapt = na,
#                n.iter= ni,
#                n.burnin = na)
# MCMCtrace(jm_red_val) 
# 
# # Process validation model output for metrics
# # Filter to validation data
# holdout <- k_df_val %>% #y_df_val to k_df_val
#   rename(full_individual=rowname) %>%
#   filter(!is.na(k_sim)) %>% 
#   rownames_to_column() %>%
#   rename(reduced_individual=rowname) %>%
#   filter(is.na(expID_val)) %>% #impID_val
#   dplyr::select(full_individual, reduced_individual, real_expID=k_sim) #, real_impID=y) #changed "predID" to "real_impID"
# 
#   # k_df_val %>% #y_df_val to k_df_val
#   # rename(full_individual=rowname) %>%
#   # filter(!is.na(k_sim)) %>% 
#   # rownames_to_column() %>%
#   # rename(reduced_individual=rowname) %>%
#   # filter(is.na(expID_val)) %>% #impID_val
#   # dplyr::select(full_individual, reduced_individual, real_expID=k_sim) #, real_impID=y) #changed "predID" to "real_impID"
# 
# # Posterior draws of predicted classifications
# full_k_out <- MCMCchains(jm_full_val, params='k') %>% #y to k
#   as_tibble() %>%
#   dplyr::select(matches(paste0("k\\[",holdout$full_individual,"\\]"))) %>% #y to k
#   mutate(draw=1:n()) %>%
#   pivot_longer(cols=-draw,names_to="full_individual",values_to="est_expID") %>% #changed from "predID" to "est_impID" to est_expID
#   mutate_at(c('full_individual'), readr::parse_number) %>%
#   mutate(full_individual=as.character(full_individual)) %>%
#   left_join(holdout %>% dplyr::select(full_individual,real_expID)) %>% #added real_impID
#   mutate(match = (est_expID==real_expID),#match = (predID==trueID),
#          model = "full")
# reduced_k_out <- MCMCchains(jm_red_val, params='k') %>%
#   as_tibble() %>%
#   #dplyr::select(matches(paste0("k\\[",holdout$reduced_individual,"\\]"))) %>%
#   mutate(draw=1:n()) %>%
#   pivot_longer(cols=-draw,names_to="reduced_individual",values_to="est_expID") %>%
#   mutate_at(c('reduced_individual'), readr::parse_number) %>%
#   mutate(reduced_individual=as.character(reduced_individual)) %>%
#   #left_join(holdout %>% dplyr::select(reduced_individual,real_expID)) %>% #added real_impID
#   mutate(match = (est_expID==real_expID),
#          model = "reduced")
# 
# Generate confusion matrices for each posterior draw ---------------------
# #max_draw <- 2000 #what does this do?
# all_combos <- expand.grid(draw = 1:max(full_k_out$draw),
#                           est_impID = 1:max(k_df_val$y), #ais changed predID_idx to predID, same from trueID
#                           trueID = 1:max(k_df_val$y))  #%>% as_tibble %>% filter(draw < max_draw) 
# 
# full_cm <- full_k_out %>%
#   #filter(draw < max_draw) %>%
#   count(draw, est_expID, real_expID) %>% #predID to est_impID
#   full_join(all_combos) %>% # fill in implicit zeros
#   reshape2::acast(draw ~ real_expID ~ est_expID, #predID to est_impID
#                   value.var = "n",
#                   fill = 0)
# reduced_cm <- reduced_k_out %>%
#   #filter(draw < max_draw) %>%
#   count(draw, est_expID, real_expID) %>%
#   full_join(all_combos) %>% # fill in implicit zeros
#   reshape2::acast(draw ~ real_expID ~ est_expID,
#                   value.var = "n",
#                   fill = 0)
# 
# assertthat::assert_that(!any(is.na(full_cm)))                # no NA vals
# assertthat::assert_that(dim(full_cm)[2] == dim(full_cm)[3])  # square matrices
# 
# get_metrics <- function(confusion_matrix) {
#   # confusion_matrix is a (true, pred) square matrix
#   true_positives <- diag(confusion_matrix)
#   false_positives <- colSums(confusion_matrix) - diag(confusion_matrix)
#   false_negatives <- rowSums(confusion_matrix) - diag(confusion_matrix)
#   precision <- true_positives / (true_positives + false_positives)
#   recall <- true_positives / (true_positives + false_negatives)
#   f1 <- 2 * (precision * recall) / (precision + recall)
#   tibble(ID = 1:length(f1),
#          precision = precision,
#          recall = recall,
#          f1 = f1)
# }
# 
# full_metrics <- apply(full_cm, 1, get_metrics) %>%
#   bind_rows(.id = "draw") %>%
#   mutate(model = "full")
# reduced_metrics <- apply(reduced_cm, 1, get_metrics) %>%
#   bind_rows(.id = "draw") %>%
#   mutate(model = "reduced") 
# 
# # Create code to list accuracy and macro-averages in one line
# val_metrics <- rbind(full_metrics, reduced_metrics) %>%
#   group_by(model,ID) %>%
#   #AIS one thing to look into for why the validaiton metrics look whack is whether I should be taking the average
#   #of a single species' metrics across draws before calculating the macro-average for the model
#   summarize(precision=mean(precision, na.rm=T), recall=mean(recall, na.rm=T), f1=mean(f1, na.rm=T)) %>%
#   #filter(ID != 3) %>%
#   group_by(model) %>%
#   summarize(precision=mean(precision, na.rm=T), recall=mean(recall, na.rm=T), f1=mean(f1, na.rm=T)) %>%
#   left_join(full_join(full_k_out,reduced_k_out) %>%
#               group_by(model) %>%
#               summarize(accuracy=mean(match))) %>%
#   mutate(fraction = fracs[current],
#          dataset = iter) %>%
#   relocate(accuracy, .after=model) %>%
#   relocate(model, .after=dataset)
# 
# val_met_iter <- rbind(val_met_iter, val_metrics)



# Goes under setting k_sim
# # Prepare input data for calculating validation metrics
# # To validate models, withhold imperfect IDs from individuals that have been validated
# #perc_withld <- 0.2 #percent individuals that have been validated that will have imperfect ID withheld
# num_withld <- round(0.1*Ltot) #take 10% of the total dataset to have imperfect IDs set to NA
# # Select rows at random that will have paratax ID withheld
# rows_withld <- data.frame(k_sim) %>%
#   rownames_to_column() %>%
#   filter(!is.na(k_sim)) %>%
#   slice_sample(n=num_withld) %>% #prop=perc_withld) %>%
#   pull(rowname) #select vector of row indices that will be withheld
# 
# # Assign parataxID_idx in selected rows to NA
# k_df_val <- data.frame(y,k_sim) %>% #y_df_val to k_df_val
#   rownames_to_column() %>%
#   mutate(expID_val = ifelse(rowname %in% rows_withld, NA, k_sim)) %>% #y to k_sim
#   # Create variable to initialize k in reduced model fit to validation dataset
#   # k.init should be NA where expID_val has a value and a randomly generated, but possible, species ID where expID_val is NA
#   mutate(k.init = ifelse(is.na(expID_val), sample(1:K_val,length(k_sim),replace=T), NA)) 

