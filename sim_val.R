# Here we generate validation metrics for our simulations 
# In our simulations, we're assessing whether



# Generate data for model input

nsite <- 30 #number of sites 
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

n <- 10 #number of fractions
fracs <- sort(runif(10, min=0.01, max=0.99), decreasing=T) #seq(0.01, 0.99, length.out = n)

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
    
    # Prepare input data for calculating validation metrics
    # To validate models, withhold imperfect IDs from individuals that have been validated
    perc_withld <- 0.2 #percent individuals that have been validated that will have imperfect ID withheld
    # Select rows at random that will have paratax ID withheld
    rows_withld <- data.frame(k_sim,y) %>%
        rownames_to_column() %>%
        filter(!is.na(k_sim)) %>%
        slice_sample(prop=perc_withld) %>%
        pull(rowname) #select vector of row indices that will be withheld
    # Require that at least one of each imperfect ID is sliced
    if (length(unique(data.frame(k_sim, y)$y)) > length(unique(data.frame(k_sim,y)$y[as.numeric(rows_withld)])) ) {
        # Swap out one of the rows_withheld indices with the imperfect ID that isn't covered
        
    }
    # Assign parataxID_idx in selected rows to NA
    y_df_val <- data.frame(k_sim, y) %>%
        rownames_to_column() %>%
        mutate(impID_val = ifelse(rowname %in% rows_withld, NA, y)) 
    
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
    na <- 6000
    
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
                        y = y_df_val$impID_val, #for all individuals, we get paratxonomist IDs
                        z = z.dat,
                        R = R)
    jm_full <- jags(data = jags_d_full,
                    inits = function(){ list(z = z.init)},
                    parameters.to.save = c("Theta", "y"), #we only need these two outputs for validation
                    model.file = "sim_full_JAGS.txt",
                    n.chains = nc,
                    n.adapt = na, 
                    n.iter= ni,
                    n.burnin = na,
                    n.thin = ni/1000)
    #saveRDS(jm_full, "output/val_full_jm.rds")
    
    # Run reduced model
    jags_d_red <- list(K_val = dim(alpha)[1],
                       K_imp = dim(alpha)[2],
                       alpha = alpha,
                       k = y_df_val %>% 
                           filter(!is.na(k_sim)) %>%
                           pull(k_sim), #true ID is known if  individual was labeled by expert
                       y = y_df_val %>%
                           filter(!is.na(k_sim)) %>%
                           pull(impID_val), #for all individuals, we get paratxonomist IDs
                       Ltot = nrow(y_df_val %>% filter(!is.na(k_sim)))) 
    jm_red <- jags(data = jags_d_red,
                   inits = NULL, 
                   parameters.to.save = c("Theta", "y"), #we only need these two outputs for validation
                   model.file = "sim_red_JAGS.txt",
                   n.chains = nc,
                   n.adapt = na, 
                   n.iter= ni,
                   n.burnin = na,
                   n.thin = ni/1000)
    #saveRDS(jm_reduced, "output/val_reduced_jm.rds")

    
    
    
    
    
    
    
    
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
    
}

# Save true params and estimates together locally 
saveRDS(jmsumm_iter, paste0("output/simulations/sim_",iter,".rds"))

    
    
    
    
    
    
    
    
    
    
    
    
    

# View JAGS output --------------------------------------------------------

full_val_jm <- readRDS("output/val_full_jm.rds")
reduced_val_jm <- readRDS("output/val_reduced_jm.rds")

### Compare the validation results

# Filter to validation data
holdout <- y_df_val %>% 
    rename(full_individual=rowname) %>%
    filter(!is.na(k_sim)) %>%
    rownames_to_column() %>%
    rename(reduced_individual=rowname) %>%
    filter(is.na(impID_val)) %>%
    dplyr::select(full_individual, reduced_individual, trueID=k_sim, predID=y)

# Posterior draws of predicted classifications
full_y_out <- MCMCchains(full_val_jm, params = 'y') %>%
    as_tibble() %>%
    dplyr::select(matches(paste0("y\\[",holdout$full_individual,"\\]"))) %>%
    mutate(draw=1:n()) %>%
    pivot_longer(cols=-draw,names_to="full_individual",values_to="predID") %>%
    mutate_at(c('full_individual'), readr::parse_number) %>%
    mutate(full_individual=as.character(full_individual)) %>%
    left_join(holdout %>% dplyr::select(full_individual, trueID)) %>%
    mutate(match = (predID==trueID), 
           model = "full")

reduced_y_out <- MCMCchains(reduced_val_jm, params = 'y') %>%
    as_tibble() %>%
    dplyr::select(matches(paste0("y\\[",holdout$reduced_individual,"\\]"))) %>%
    mutate(draw=1:n()) %>%
    pivot_longer(cols=-draw,names_to="reduced_individual",values_to="predID") %>%
    mutate_at(c('reduced_individual'), readr::parse_number) %>%
    mutate(reduced_individual=as.character(reduced_individual)) %>%
    left_join(holdout %>% dplyr::select(reduced_individual, trueID)) %>%
    mutate(match = (predID==trueID), 
           model = "reduced")


# Generate confusion matrices for each posterior draw ---------------------
max_draw <- 2000
all_combos <- expand.grid(draw = 1:max(full_y_out$draw), 
                          predID_idx = 1:max(y_df_val$y), 
                          trueID_idx = 1:max(y_df_val$y))  %>%
    as_tibble %>%
    filter(draw < max_draw)

full_cm <- full_y_out %>%
    filter(draw < max_draw) %>%
    count(draw, predID, trueID) %>%
    full_join(all_combos) %>% # fill in implicit zeros
    reshape2::acast(draw ~ trueID ~ predID, 
                    value.var = "n", 
                    fill = 0)

reduced_cm <- reduced_y_out %>%
    filter(draw < max_draw) %>%
    count(draw, predID, trueID) %>%
    full_join(all_combos) %>% # fill in implicit zeros
    reshape2::acast(draw ~ trueID ~ predID, 
                    value.var = "n", 
                    fill = 0)

assertthat::assert_that(!any(is.na(full_cm)))                # no NA vals
assertthat::assert_that(dim(full_cm)[2] == dim(full_cm)[3])  # square matrices
if (dim(full_cm)[2] != dim(full_cm)[3]) {
    full_cm_temp <- array(NA, dim = c(dim(full_cm)[1],length(unique(y_df_val$y)),length(unique(y_df_val$y)))) #initialize
    reduced_cm_temp <- array(NA, dim = c(dim(reduced_cm)[1],length(unique(y_df_val$y)),length(unique(y_df_val$y)))) #initialize
    for (i in 1:dim(full_cm)[3]) {
        # Add a column for the 3rd imperfect ID
        #full_cm_temp[,,i] <- cbind(full_cm[,,i], rep(0,dim(full_cm)[1]))
        reduced_cm_temp[,,i] <- cbind(reduced_cm[,,i], rep(0,dim(reduced_cm)[1]))
    }
    full_cm <- full_cm_temp
    reduced_cm <- reduced_cm_temp
}
assertthat::assert_that(dim(full_cm)[2] == dim(full_cm)[3])  # square matrices

get_metrics <- function(confusion_matrix) {
    # confusion_matrix is a (true, pred) square matrix
    true_positives <- diag(confusion_matrix)
    false_positives <- colSums(confusion_matrix) - diag(confusion_matrix)
    false_negatives <- rowSums(confusion_matrix) - diag(confusion_matrix)
    precision <- true_positives / (true_positives + false_positives)
    recall <- true_positives / (true_positives + false_negatives)
    f1 <- 2 * (precision * recall) / (precision + recall)
    tibble(impID = 1:length(f1),
           precision = precision, 
           recall = recall, 
           f1 = f1)
}

full_metrics <- apply(full_cm, 1, get_metrics) %>%
    bind_rows(.id = "draw") %>%
    mutate(model = "full")
reduced_metrics <- apply(reduced_cm, 1, get_metrics) %>%
    bind_rows(.id = "draw") %>%
    mutate(model = "reduced")

# Holdout log likelihood (log probabilities from the categorical distribution) 
# (could use dcat() function: https://www.rdocumentation.org/packages/mvc/versions/1.3/topics/dcat)
theta_summ_full_val <- MCMCchains(full_val_jm, params = 'Theta') %>% 
    as_tibble() %>%
    mutate(draw=1:n()) %>%
    pivot_longer(cols=-draw,names_to="index",values_to="value") %>%
    separate("index", into = c("trueID", "predID"), sep = ",") %>%
    mutate_at(c('trueID', 'predID'), readr::parse_number) 
theta_summ_reduced_val <- MCMCchains(reduced_val_jm, params = 'Theta') %>% 
    as_tibble() %>%
    mutate(draw=1:n()) %>%
    pivot_longer(cols=-draw,names_to="index",values_to="value") %>%
    separate("index", into = c("trueID", "predID"), sep = ",") %>%
    mutate_at(c('trueID', 'predID'), readr::parse_number) 

full_loglik <- holdout %>%
    count(trueID, predID) %>%
    left_join(theta_summ_full_val) %>%
    group_by(draw) %>%
    summarize(log_lik = sum(n * log(value))) %>%
    mutate(model = "full")
reduced_loglik <- holdout %>%
    count(trueID, predID) %>%
    left_join(theta_summ_reduced_val) %>%
    group_by(draw) %>%
    summarize(log_lik = sum(n * log(value))) %>%
    mutate(model = "reduced")


# create validation.png: accuracy, f1 score, holdout log-lik
acc <- full_y_out %>% 
    group_by(draw) %>%
    summarize(accuracy=mean(match)) %>%
    mutate(model="full") %>%
    rbind(reduced_y_out %>% 
              group_by(draw) %>%
              summarize(accuracy=mean(match)) %>%
              mutate(model="reduced")) %>%
    ggplot(aes(x=accuracy, fill=model)) +
    geom_density(alpha = .6) + 
    xlab("Accuracy") +
    ylab("Density") +
    theme(legend.position = "none")     

f1 <- full_join(full_metrics, reduced_metrics) %>%
    group_by(draw, model) %>%
    summarize(macro_f1 = mean(f1, na.rm = TRUE)) %>%
    ggplot(aes(macro_f1, fill = model)) + 
    geom_density(alpha = .6) +
    xlab("F1 score") +
    ylab(NULL) + 
    theme(legend.position = "none")     

hll <- full_join(full_loglik, reduced_loglik) %>%
    ggplot(aes(log_lik, fill = model)) + 
    geom_density(alpha = .6) + 
    xlab("Holdout log-likelihood") + 
    ylab(NULL) + 
    theme(legend.position = "none")   

legend <- get_legend( full_join(full_loglik, reduced_loglik) %>%
                          ggplot(aes(log_lik, fill = model)) + 
                          geom_density(alpha=0.6) +
                          scale_fill_discrete(labels = c("full model", "reduced model")) +
                          theme(legend.key.size = unit(0.8, 'cm'),
                                legend.title = element_blank(),
                                legend.text = element_text(size=10)) )

grid.arrange(acc, f1, hll, legend, ncol = 4) 
val <- arrangeGrob(acc, f1, hll, legend, ncol = 4)
ggsave("figures/validation.png", val, width = 8)
