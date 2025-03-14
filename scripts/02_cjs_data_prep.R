# Purpose: Prepare capture history data for CJS survival modeling with program
# Mark and Jags
# Author: Ryan N. Kinzer
# Date: 2021-10-08

# load pkgs----
library(tidyverse)

# load data----
sy <- 2024
load(paste0('./data/kelt_data_sy',sy,'.Rda'))

# create model data----
ch_mod_dat <- tag_summary %>%
  select(-other) %>% # extra field from previous step
  filter(spawner_above == 1) %>%
  #filter(!(POP_NAME %in% c('GRA', 'RAPH', 'OXBO', 'DNFH'))) %>% #hatcheries and unknown spawning areas
  rowwise() %>%
  mutate(kelt_below = max(c_across(c(kelt_GOA:rs_above))),
         kelt_rs = max(c_across(c(kelt_BON, rs_BON, rs_LGR, rs_above))),
         rs_all = max(c_across(c(rs_BON, rs_LGR, rs_above)))) %>%
  mutate(rs_LGR = max(c_across(c(rs_LGR, rs_above)))) %>%
  ungroup()

tmp <- ch_mod_dat %>%
  group_by(spawn_yr) %>%
  summarise(across(.cols = c(spawner_above:rs_all, ), sum))
  
save(ch_mod_dat, file = paste0('./data/cjs_model_data_sy',sy,'.Rda'))

# select parameters
# pops <- c('GRJOS-s', 'IRMAI-s', 'GRUMA-s')
# yr <- 2011:2020
# 
# # subset data
# ch_filtered <- ch_mod_dat %>%
#   ungroup() %>%
#   filter(TRT_POPID %in% pops) %>%
#   filter(spawn_year %in% yr)
# 
# ch_filtered %>%
#   group_by(TRT_POPID, spawn_year) %>%
#   summarise(across(.cols = c(spawner_above:rs_above), sum))
