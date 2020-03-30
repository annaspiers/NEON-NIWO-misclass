# This script builds a dynamic occupancy model for a single species using simulated data. A dynamic occupancy model estimates changes in z, probability of occupancy, across seasons. Followed model in section 13.5 of Bayesian Population Analysis (Kery & Schaub 2012).

# Seasons over a series of years represent an common, natural grouping factor, where each year is the primary sampling occasion and the surveys within each year are secondary sampling occasions. It is natural then to assume closure among secondary seasons within a primary season.


# Generate data -------------
# Function to simulate detection/nondetection data for dynamic site-occ model. Annual variation in probabilities of patch survival, colonization and detection is specified by the bounds of a uniform distribution. 
data.fn <- function( R = 250,  #number of sites across a population
                     J = 3, #number of replicate surveys
                     K= 10, #number of survey seasons
                     psi1 = 0.4, #occupancy probability in first year 
                     range.p = c(0.2, 0.4), #bounds of uniform dist from which annual p is drawn 
                     range.phi = c(0.6, 0.8), #bounds of uniform dist from which survival prob is drawn 
                     range.gamma = c(0, 0.1) #bounds of uniform dist from which colonization prob is drawn 
                     ) { 

    # Set up some required arrays 
    site <- 1:R 
    year <- 1:K 
    psi <- rep(NA, K) # Occupancy probability
    muZ <- z <- array(dim = c(R, K)) # Expected and realized occurrence 
    y <- array(NA, dim = c(R, J, K)) # Detection histories
    
    # Determine initial occupancy and demographic parameters 
    psi[1] <- psi1 # Initial occupancy probability
    p <- runif(n = K, min = range.p[1], max = range.p[2]) 
    phi <- runif(n = K-1, min = range.phi[1], max = range.phi[2]) 
    gamma <- runif(n = K-1, min = range.gamma[1], max = range.gamma[2])
    
    # Generate latent states of occurrence 
    # First year 
    z[,1] <- rbinom(R, 1, psi[1]) # Initial occupancy state
    
    # Later years 
    for(i in 1:R){ # Loop over sites
        for(k in 2:K){ # Loop over years
            muZ[k] <- z[i, k-1]*phi[k-1] + (1-z[i, k-1])*gamma[k-1] # Prob for occ.
            z[i,k] <- rbinom(1, 1, muZ[k]) 
        }#k
    }#i
    
    # Plot realised occupancy 
    plot(year, apply(z, 2, mean), type = "l", xlab = "Year", ylab = "Occupancy or Detection prob.", 
         col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1) 
    lines(year, p, type = "l", col = "red", lwd = 2, lty = 2)
    
    # Generate detection/nondetection data 
    for(i in 1:R){ 
        for(k in 1:K){ 
            prob <- z[i,k] * p[k] 
            for(j in 1:J){ 
                y[i,j,k] <- rbinom(1, 1, prob) 
            }#j
        }#k
    }#i
    
    # Compute annual population occupancy 
    for (k in 2:K){ 
        psi[k] <- psi[k-1]*phi[k-1] + (1-psi[k-1])*gamma[k-1] 
    }
    
    # Plot apparent occupancy 
    psi.app <- apply(apply(y, c(1,3), max), 2, mean) 
    lines(year, psi.app, type = "l", col = "black", lwd = 2) 
    text(0.85*K, 0.06, labels = "red solid – true occupancy\n red dashed – detection probability\n black – observed occupancy")
    
    # Return data 
    return(list(R = R, J = J, K= K, psi = psi, psi.app = psi.app, z = z, phi = phi, gamma = gamma, p = p, y = y))
}

data <- data.fn() #AIS added in these two lines to make variables from data.fn() readable - I wonder how they did it without
list2env(data, envir = .GlobalEnv)

# Bundle data
str(JAGSdata <- list(y = y, 
                     nsite = dim(y)[1], 
                     nrep = dim(y)[2], 
                     nyear = dim (y)[3])) 
# Initial values 
zst <- apply(y, c(1, 3), max) 
# Observed occurrence as inits for z 
JAGSinits <- function(){ list(z = zst)} 
# Parameters monitored 
JAGSparams <- c("psi", "phi", "gamma", "p", "n.occ", "growthr", "turnover")

# MCMC settings
nc <- 3 #MCMC chains         
ni <- 2500 #MCMC iterations
nb <- 500 #MCMC burnin
nt <- 4  #MCMC thin

# JAGS model
JAGSout <- jags(data = JAGSdata,
                 inits = JAGSinits,
                 parameters.to.save = JAGSparams,
                 model.file = "occupancy/03_bpa_dynamic_single-sp_JAGS.txt",
                 n.chains = nc,
                 n.iter = ni,
                 n.burnin = nb,
                 n.thin = nt)

print(JAGSout, dig=2) #growthr <1 means the number of occupied sites is smaller between years
names(JAGSout)
plot(JAGSout) 

# Compare data value to credible interval range
# Occupancy probability
print(cbind(data$psi, JAGSout$summary[1:K, c(1, 2, 3, 7)]), dig = 3) 
# Survival probability
print(cbind(data$phi, JAGSout$summary[(K+1):(K+(K-1)), c(1, 2, 3, 7)]), dig = 3) 
# Colonization probability
print(cbind(data$gamma, JAGSout$summary[(2*K):(2*K+(K-2)), c(1, 2, 3, 7)]), dig = 3) 
# Detection probability
print(cbind(data$p, JAGSout$summary[(3*K-1):(4*K-2), c(1, 2, 3, 7)]), dig = 3)

# Plot occupancy through seasons
plot(1:K, data$psi, type = "l", xlab = "Year", ylab = "Occupancy probability", col = "red", xlim = c(0,K+1), ylim = c(0,1), lwd = 2, lty = 1, frame.plot = FALSE, las = 1)
lines(1:K, data$psi.app, type = "l", col = "black", lwd = 2) 
points(1:K, JAGSout$mean$psi, type = "l", col = "blue", lwd = 2) 
segments(1:K, JAGSout$summary[1:K,3], 1:K,JAGSout$summary[1:K,7], col = "blue", lwd = 1)
