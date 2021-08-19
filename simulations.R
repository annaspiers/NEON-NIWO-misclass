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


# Simulate data from model priors
# By generating data from the model, fitting the model to these data, then computing posteriors for (occupancy, demographic rate) parameter estimates, we have theoretical guarantees about interval coverage. 

# Generate data from model priors --------------------------------------------------

#iter <- 1
#fracs <- c(0.8, 0.3)

sim_1dataset <- function(iter, fracs) {
  
  nsite <- 20 #number of sites #bump back up to 30 for final processing
  nsurv <- 3  #number of suveys in the season
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
  diag(alpha) <- 500 # couldn't find great match with 2/3 species here to 
  #diag(alpha) <- 75 #AIS used to be 2/200, 2/75 and 2/40 had similar results
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
  
  # Pi (dirichlet distribution for each encountered individual)
  # k (validated species for an individual)
  # y (imperfectly observed state for an individual)
  Ltot <- sum(L)
  site <- sample(1:nsite, Ltot, replace=TRUE)             #AIS is this workaround ok for initializing site?
  pi <- matrix(NA, nrow=Ltot, ncol=K_val)
  k <- y <- rep(NA, Ltot)
  for (l in 1:Ltot) {
    pi[l, 1:K_val] <- zlam[site[l], 1:K_val] / sum(zlam[site[l], 1:K_val])
    #if (is.na(pi[l,1])) { #arbitrary column since all in a row would be NA
    #pi[l, ] <- rep(0,K_val)
    
    #pi[l, ] <- rep(0,K_val) #assign all values to 0
    #random_ind <- sample(K_val,1)
    #pi[l, random_ind] <- 1 #but assign one of them 1 randomly
    
    #pi[l, ] <- rep(1/K_val,K_val)
    #} #AIS is this an ok workaround? When pi has a row of NA's, then k is NA for that row, which we don't want. Instead, I give each K_val the same weight
    k[l] <- rcat(1, pi[l, 1:K_val])
    y[l] <- rcat(1, Theta[k[l], 1:K_imp])
  }
  
  
  # Vary the fraction of individuals that are validated ---------------------
  
  # For this simulation, we are changing the fraction of samples that were validated.
  # At this point, the object k stores the validated state value for 100% of samples.
  # We need to vary the the number of values in k that are NA, thereby varying the fraction
  # of individuals that are validated 
  
  k_mat <- matrix(data=k, nrow=Ltot, ncol=length(fracs))
  colnames(k_mat) <- paste0("k",fracs)
  for (j in 1:ncol(k_mat)) {
    if (j!=1) {k_mat[,j] <- k_mat[,j-1]}
    alreadyNA <- sum(is.na(k_mat[,j]))
    needtobeNA <- round((1-fracs[j])*Ltot)
    ind <- sample(which(!is.na(k_mat[,j])), (needtobeNA - alreadyNA), replace=F)
    k_mat[ind,j]<- NA
  }
  
  # Initialize df to store jmsummm for entire dataset (includes all k values)
  jmsumm_iter <- data.frame()
  
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
    z.dat <- array(NA, dim = c(nsite, K_val),
                   dimnames = list(1:nsite, 1:K_val))
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
    ni <- 30000
    na <- 5000
    
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
                    inits = function(){ list(z = z.init)},
                    parameters.to.save = c("psi","lambda", "Theta"), 
                    model.file = "sim_full_JAGS.txt",
                    n.chains = nc,
                    n.adapt = na, 
                    n.iter= ni,
                    n.burnin = na,
                    n.thin = ni/1000)
    # AIS temporary while troubleshooting
    #saveRDS(jm_full, paste0("output/simulations/sim_3_jm_full_0.001.rds"))
    #saveRDS(jm_full, paste0("output/simulations/sim_3_jm_full_0.05.rds"))
    #saveRDS(jm_full, paste0("output/simulations/sim_3_jm_full_1.0.rds"))
    # now compare posterior
    
    jmsumm_full <- MCMCsummary(jm_full) %>%
      mutate(fraction = fracs[current],
             dataset = iter,
             model = "full") %>%
      rownames_to_column()
    
    
    
    # Run reduced model
    jags_d_red <- list(K_val = dim(alpha)[1],
                       K_imp = dim(alpha)[2],
                       alpha = alpha,
                       k = data.frame(y= y, k = k_mat[ ,current]) %>% 
                         filter(!is.na(k)) %>%
                         pull(k), #true ID is known if  individual was labeled by expert
                       y = data.frame(y= y, k = k_mat[ ,current]) %>%
                         filter(!is.na(k)) %>%
                         pull(y), #for all individuals, we get paratxonomist IDs
                       Ltot = nrow(data.frame(y= y, k = k_mat[ ,current]) %>%
                                     filter(!is.na(k)))) 
    jm_red <- jags(data = jags_d_red,
                   inits = NULL, 
                   parameters.to.save = c("Theta"),
                   model.file = "sim_red_JAGS.txt",
                   n.chains = nc,
                   n.adapt = na, 
                   n.iter= ni,
                   n.burnin = na,
                   n.thin = ni/1000)
    jmsumm_red <- MCMCsummary(jm_red) %>%
      mutate(fraction = fracs[current],
             dataset = iter,
             model = "reduced") %>%
      rownames_to_column()
    
    # Iniitalize dataframes that store true parameter values
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
    
    # Write cluster status to a txt file
    write(paste0("Finished dataset ",iter," with ",fracs[current]," samples validated. Took ",(Sys.time() - start_time)," (",Sys.time(),")"),file="monitoring.txt",append=TRUE)
  }
  
  # Save true params and estimates together locally 
  saveRDS(jmsumm_iter, paste0("output/simulations/sim_",iter,".rds"))
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
iter_tot <- 2
system.time( 
  clusterApply(cl, 
               1:iter_tot, #iter
               sim_1dataset, 
               c(0.8, 0.2)) #c(1.0, 0.9, 0.75, 0.5, 0.25, 0.15))  #fracs
)
# use top to monitor progress

stopCluster(cl)




### SCRAP NOTES

# Also, save your results now and again so you at least have something if your
# program crashes after many hours of work.
# save(myresults,file="results.Rdata")
#start_time <- Sys.time(); end_time <- Sys.time(); total_time <- end_time - start_time

