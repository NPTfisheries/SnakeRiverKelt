# Purpose: Basic CJS model for steelhead kelt survival.
# Author: Ryan N. Kinzer
# Date: 2021-10-08

cjsBasic <- "model {    
  for (i in 1:nind){                       #loops through each individual
    for (t in f[i]:(n.occasions-1)){      #loops through each capture period
      phi[i,t] <- Sur[t]                 #relates trap survival to individual
      p[i,t] <- Det[t]                   #relates trap detection to individual
    } #t
  } #i
  
  for (t in 1:(n.occasions-1)){
    Sur[t] ~ dunif(0, 1)        # Priors for time-spec. survival
    Det[t] ~ dunif(0, 1)         # Priors for time-spec. detection
  }
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
  } #i
}"

cat(cjsBasic, file = './modelfiles/cjsBasic.txt')
