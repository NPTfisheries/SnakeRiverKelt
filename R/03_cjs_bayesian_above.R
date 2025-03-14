# Prepare data and run Bayesian CJS models.
# Author: Ryan N. Kinzer
# Date: 2021-10-08

#library(runjags)
#library(rjags)
library(tidyverse)
library(jagsUI)
source('./R/run_cjs_bayes.R')

# load data----

yr <- 2024

load(paste0('./data/cjs_model_data_',yr,'.Rda'))

mod_dat <- ch_mod_dat %>%
  filter(spawner_above == 1) %>%
  filter(!is.na(TRT_POPID)) %>%
  select(spawner_above, kelt_above, kelt_LGR, kelt_rs, group1 = spawn_year, group2 = TRT_POPID)

# run JAGS models----
iters <- 20000
burnin <- 5000
chains <- 3
thin <- 3
output_file <- paste0('./data/mod_runs_above_',yr,'.Rda')

run_cjs_bayes(data = mod_dat, iters = iters, burnin = burnin, chains = chains, thin = thin, output_file = output_file)

