# Purpose: Define CJS models for kelt phi.tvival estimates
# Author: Ryan N. Kinzer
# Created: 2021-10-14

# Cormack-Jolly-Seber Model----
# assumes all spawn years and populations are equal with fixed time effects
cjsFixedTime <- 'model{    
  
  # priors
  for (t in 1:(n.occasions-1)){
    PHI[t] ~ dunif(0, 1)        # time specific survival
    P[t] ~ dunif(0, 1)         # time specific detection
  }
  
  # likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(mu1[i,t])
      mu1[i,t] <- PHI[t-1] * z[i,t-1]
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- P[t-1] * z[i,t]
    } #t
  } #i
}'

cat(cjsFixedTime, file = './modelfiles/cjsFixedTime.txt')


# Random effect group and fixed time effects Cormack-Jolly-Seber Model----
# assumes spawn year survivals are different -- good for a single population.
cjsFixedTimeRandomGroup <- 'model{    
  for (i in 1:nind){               
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- PHI[group1[i], t]  
      p[i,t] <- P[group1[i], t]      
    } #t
  } #i
  
  for (t in 1:(n.occasions-1)){
    for(l in 1:g1){
      logit(PHI[l, t]) <- mu.phi[t] + group1.phi[l, t]
      logit(P[l, t]) <- mu.p[t] + group1.p[l, t]
      group1.phi[l, t] ~ dnorm(0,tau.group1.phi)
      group1.p[l, t] ~ dnorm(0,tau.group1.p)         
    } #l  
  } #t
  
  for(t in 2:n.occasions-1){
    mu.phi[t] ~ dnorm(0,.1)        
    mu.p[t] ~ dnorm(0,.1)
    mean.phi[t] <- exp(mu.phi[t])/(1+exp(mu.phi[t]))
    mean.p[t] <- exp(mu.p[t])/(1+exp(mu.p[t]))
  }
  
  sigma.group1.phi ~ dunif(0,1)            
  tau.group1.phi <- 1/(sigma.group1.phi*sigma.group1.phi)
  sigma.group1.p ~ dunif(0,1)             
  tau.group1.p <- 1/(sigma.group1.p*sigma.group1.p)
  
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
}'

cat(cjsFixedTimeRandomGroup, file = './modelfiles/cjsFixedTimeRandomGroup.txt')

# Two random group effect Cormack-Jolly-Seber Model----


cjsFixedTimeRandomGroups <- 'model{    
  # Constraints
  for (i in 1:nind){               
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- PHI[group2[i], group1[i], t]  
      p[i,t] <- P[group2[i], group1[i], t]      
    } #t
  } #i
  
  for (t in 1:(n.occasions-1)){
    for(l in 1:g1){
      for(k in 1:g2){
        logit(PHI[k,l,t]) <- mu.phi[t] + group1.phi[l,t] + group2.phi[k,t]
        logit(P[k,l,t]) <- mu.p[t] + group1.p[l,t] + group2.p[k,t]
      } #k
    } #l  
  } #t
  
  # Priors for mean time effects
  for(t in 2:n.occasions-1){
    mu.phi[t] ~ dnorm(0,.1)        
    mu.p[t] ~ dnorm(0,.1)
    mean.phi[t] <- exp(mu.phi[t])/(1+exp(mu.phi[t]))
    mean.p[t] <- exp(mu.p[t])/(1+exp(mu.p[t]))
  }
  
  # Priors for random group 1 effects
  for (t in 1:(n.occasions-1)){
   for(l in 1:g1){ 
    group1.phi[l,t] ~ dnorm(0,tau.group1.phi)
    group1.p[l,t] ~ dnorm(0,tau.group1.p)
   } #l  
  } #t
  
  sigma.group1.phi ~ dunif(0,1)            
  tau.group1.phi <- 1/(sigma.group1.phi*sigma.group1.phi)
  sigma.group1.p ~ dunif(0,1)             
  tau.group1.p <- 1/(sigma.group1.p*sigma.group1.p)
  
  # Priors for random group 2 effects
  for (t in 1:(n.occasions-1)){
   for(k in 1:g2){ 
    group2.phi[k,t] ~ dnorm(0,tau.group2.phi)
    group2.p[k,t] ~ dnorm(0,tau.group2.p)
   } #k  
  } #t
  
  sigma.group2.phi ~ dunif(0,1)            
  tau.group2.phi <- 1/(sigma.group2.phi*sigma.group2.phi)
  sigma.group2.p ~ dunif(0,1)             
  tau.group2.p <- 1/(sigma.group2.p*sigma.group2.p)
  
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
}'

cat(cjsFixedTimeRandomGroups, file = './modelfiles/cjsFixedTimeRandomGroups.txt')