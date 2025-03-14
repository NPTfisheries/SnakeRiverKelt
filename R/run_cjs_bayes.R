# Purpose: function to run three bayesian CJS models
# Author: Ryan N. Kinzer
# Created: 2021
# Modified: 31 May 2024

run_cjs_bayes <- function(data, iters, burnin, chains, thin, output_file){
  
  start <- Sys.time()
  
  source('./R/cjs_prep_functions.R')
  source('./R/cjs_models.R')
  
  # model data
  CH <- as.matrix(data[, sapply(data, class) == "numeric"])
  f <- apply(CH, 1, get.first)
  z <- known.state.cjs(CH)
  
  # time difference model
  cjs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], 
                   z = known.state.cjs(CH))
  
  # inits
  inits<-function(){list(z=cjs.init.z(CH,f),PHI=runif(dim(CH)[2]-1,0,1),P=runif(dim(CH)[2]-1,0,1))}
  initial<-list(inits(),inits(),inits())
  
  # set parameters to save
  jagsParams = c('PHI', 'P')
  
  mod.t = jagsUI(data = cjs.data,
                 parameters.to.save = jagsParams,
                 model.file = './modelfiles/cjsFixedTime.txt',
                 n.chains = chains,
                 n.iter = iters,
                 n.burnin = burnin,
                 n.thin = thin,
                 verbose = F,
                 parallel=TRUE)
  
  # random year----
  grp1 <- levels(as.factor(data$group1))
  group1 <- as.integer(as.factor(data$group1))
  g1 = length(unique(group1))
  
  cjs.data <- list(y = CH, group1 = group1, g1 = g1,
                   f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], 
                   z = known.state.cjs(CH))
  
  # inits
  inits<-function(){list(z=cjs.init.z(CH,f),
                         mu.phi = runif(dim(CH)[2]-1,-1,1),             # OR of mean survial
                         mu.p = runif(dim(CH)[2]-1,-1,1),
                         group1.phi = matrix(runif(((dim(CH)[2]-1)*length(unique(group1))),-1,1),
                                             nrow=length(unique(group1)),
                                             ncol=dim(CH)[2]-1),         # OR of yearly random survival effects 
                         group1.p = matrix(runif(((dim(CH)[2]-1)*length(unique(group1))),-1,1),
                                           nrow=length(unique(group1)),
                                           ncol=dim(CH)[2]-1),         # OR of yearly random detection effects
                         sigma.group1.phi = runif(1,0,1),                        # prior for survival variability
                         sigma.group1.p = runif(1,0,1))}                         # prior for detection variability
  
  initial<-list(inits(),inits(),inits())
  
  # parameters to monitor
  jagsParams = c("mean.phi","mean.p",
                 'PHI', 'P',
                 "sigma.group1.phi","sigma.group1.p")
  
  # run JAGS model
  mod.t.year = jagsUI(data = cjs.data,
                      parameters.to.save = jagsParams,
                      model.file = './modelfiles/cjsFixedTimeRandomGroup.txt',
                      n.chains = chains,
                      n.iter = iters,
                      n.burnin = burnin,
                      n.thin = thin,
                      verbose = F,
                      parallel = TRUE)
  
  ## random population and year----
  grp2 <- levels(as.factor(data$group2))
  group2 <- as.integer(as.factor(data$group2))
  g2 <- length(unique(group2))
  
  cjs.data <- list(y = CH, group1 = group1, g1 = g1,
                   group2 = group2, g2 = g2,
                   f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], 
                   z = known.state.cjs(CH))
  
  # inits
  inits<-function(){list(z=cjs.init.z(CH,f),
                         mu.phi = runif(dim(CH)[2]-1,-1,1),             # OR of mean survial
                         mu.p = runif(dim(CH)[2]-1,-1,1),
                         group1.phi = matrix(runif(((dim(CH)[2]-1)*length(unique(group1))),-1,1),
                                             nrow=length(unique(group1)),
                                             ncol=dim(CH)[2]-1),         # OR of yearly random survival effects 
                         group1.p = matrix(runif(((dim(CH)[2]-1)*length(unique(group1))),-1,1),
                                           nrow=length(unique(group1)),
                                           ncol=dim(CH)[2]-1),         # OR of yearly random detection effects
                         sigma.group1.phi = runif(1,0,1),                        # prior for survival variability
                         sigma.group1.p = runif(1,0,1),
                         
                         group2.phi = matrix(runif(((dim(CH)[2]-1)*length(unique(group2))),-1,1),
                                             nrow=length(unique(group2)),
                                             ncol=dim(CH)[2]-1),         # OR of yearly random survival effects 
                         group2.p = matrix(runif(((dim(CH)[2]-1)*length(unique(group2))),-1,1),
                                           nrow=length(unique(group2)),
                                           ncol=dim(CH)[2]-1),         # OR of yearly random detection effects
                         sigma.group2.phi = runif(1,0,1),                        # prior for survival variability
                         sigma.group2.p = runif(1,0,1))}
  
  initial<-list(inits(),inits(),inits())
  
  # parameters to monitor
  jagsParams = c('mean.phi', 'mean.p',
                 'PHI', 'P',
                 'sigma.group1.phi', 'sigma.group1.p',
                 'sigma.group2.phi', 'sigma.group2.p'
  )
  
  # run JAGS model
  mod.t.year.pop = jagsUI(data = cjs.data,
                          parameters.to.save = jagsParams,
                          model.file = './modelfiles/cjsFixedTimeRandomGroups.txt',
                          n.chains = chains,
                          n.iter = iters,
                          n.burnin = burnin,
                          n.thin = thin,
                          verbose = F,
                          parallel = TRUE)
  
  end <- Sys.time()

  
  model_ls <- list("models" = list("mod.t" = mod.t,
                                 "mod.t.year" = mod.t.year,
                                 "mod.t.year.pop" = mod.t.year.pop),
                 "data" = list("mod_dat" = data,
                               "cjs.data" = cjs.data,
                               "grp1" = grp1,
                               "grp2" = grp2),
                 "meta" = list("start" = start,
                               "end" = end,
                               "iters" = iters,
                               "burnin" = burnin,
                               "chains" = chains,
                               "thin" = thin)
                 )
  
  save(model_ls,
       file = output_file)
  
}