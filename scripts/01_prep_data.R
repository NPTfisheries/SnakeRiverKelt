# Purpose: Load data from SnakeRiverFishStatus project and develop inputs for
#   steelhead spawner, kelt, and repeat spawner summaries.
#
# Author: Ryan N. Kinzer
# Date Created: 2021-10-08
#   Date Last Modified: 2025-03-17
#   Modified By: Mike Ackerman

# load packages
library(tidyverse)
library(sf)

# source some helper functions
source('./R/find_max_spawner.R')
source("./R/steelheadCapHist.R")

# set some parameters
spp = "Steelhead"
yr_range = 2010:2024
max_sy = max(yr_range)

#---------------
# load data

# set some file paths to datasets
snake_proj = "../SnakeRiverFishStatus/"
PITcleanr_folder = paste0(snake_proj, "output/PITcleanr/human_reviewed/")
lh_folder = paste0(snake_proj, "output/life_history/")
trap_path = paste0(snake_proj, "data/LGTrappingDB/LGTrappingDB_2024-12-26.csv")

# load configuration files (just need crb_sites_sf object)
load(file = paste0(snake_proj, "data/configuration_files/site_config_LGR_20241226.rda")) ; rm(flowlines, parent_child, configuration, sr_site_pops)

# load PITcleanr cleaned steelhead observation data
pitcleanr_df = list.files(path = PITcleanr_folder,
                          pattern = "Steelhead",
                          full.names = T) %>%
  map_df(~ {
    # extract the spawn year from the file name
    spawn_yr = str_extract(.x, "(?<=SY)\\d{4}")
    readxl::read_xlsx(path = .x, sheet = "Sheet1") %>%
      mutate(spawn_yr = as.numeric(spawn_yr))
  })

# load life history data
tag_meta = list.files(path = lh_folder,
                      pattern = "Steelhead",
                      full.names = T) %>%
  map_df(~ readxl::read_xlsx(
    path = .x,
    sheet = "tag_lh"
  )) %>%
  rename(spawn_yr = spawn_year) %>%
  janitor::clean_names()

# load lgtrappingdb
# trap_df = read_csv(trap_path) %>%
#   # some LOSALM & LOCLWR genstocks are errantly coded as LOWSALM & LOWCLWR
#   mutate(GenStock = recode(GenStock, "LOWSALM" = "LOSALM")) %>%
#   mutate(GenStock = recode(GenStock, "LOWCLWR" = "LOCLWR"))

#---------------
# generate a common complete tag history dataset
cth_df = pitcleanr_df %>%
  # remove spawner observations considered FALSE for dabom
  filter(!(life_stage == "spawner" & user_keep_obs == FALSE)) %>%
  # create site code from node
  mutate(site_code = str_remove(node, "_[D|U]")) %>%
  # grab the date-time that each fish last left LGR
  group_by(spawn_yr, tag_code) %>%
  mutate(lgr_max_det = max(max_det[site_code == "LGR" & life_stage == "spawner"])) %>%
  ungroup() %>%
  # join rkm and rkm_total for each site
  left_join(crb_sites_sf %>%
              select(site_code,
                     rkm,
                     rkm_total) %>%
              st_drop_geometry(),
            by = "site_code") %>%
  # a little cleaning
  select(spawn_yr,
         id,
         tag_code,
         life_stage,
         site_code,
         node,
         rkm,
         rkm_total,
         direction,
         slot,
         min_det,
         max_det,
         lgr_max_det,
         path)

# use find_max_spawner to fix errant kelt observations?
fix_errant_kelt_calls = T

# apply find_max_spawner only if fix_errant_kelt_calls is TRUE
if (fix_errant_kelt_calls) {
  cth_df = cth_df %>%
    group_by(spawn_yr) %>%
    nest() %>%
    mutate(new = map(data, find_max_spawner)) %>%
    select(-data) %>%
    unnest(new) 
}

# RK direction
rk_ch <- cth_df %>%
  ungroup() %>%
  group_by(spawn_yr) %>%
  nest() %>%
  mutate(ch = map(data,
                  .f = ~steelheadCapHist(.))) %>%
  select(-data) %>%
  unnest(ch)

# MA direction
source("./R/steelheadCapHist_MA.R")
ma_ch = cth_df %>%
  ungroup() %>%
  group_by(spawn_yr, tag_code) %>%
  steelheadCapHist_MA() %>%
  ungroup() %>%
  select(-lgr_max_det, -grs_kelt_det)

# QA/QC: summed cap hists by spawn year
rk_ch_summary = rk_ch %>%
  group_by(spawn_yr) %>%
  summarise(across(-tag_code, sum)) %>%
  ungroup() %>%
  select(-other)

# QA/QC: summed cap hists by spawn year
ma_ch_summary = ma_ch %>%
  select(spawn_yr, 
         tag_code,
         release_lgr:rs_above) %>%
  group_by(spawn_yr) %>%
  summarise(across(-tag_code, sum)) %>%
  ungroup()

# join ma & rk capture histories
ch_compare_df = ma_ch %>%
  select(spawn_yr:rs_above) %>%
  mutate(user = "ma") %>%
  bind_rows(rk_all %>%
              select(-other) %>%
              mutate(user = "rk")) %>%
  select(spawn_yr, tag_code, user, everything()) %>%
  arrange(spawn_yr, tag_code, user)

# return records where capture histories don't match (takes awhile)
mismatch_ch = ch_compare_df %>%
  group_by(spawn_yr, tag_code) %>%
  filter(any(across(release_lgr:rs_above, ~ . != first(.)))) %>%
  ungroup()

# quick checks  
nrow(mismatch_ch) / 2                           # number of fish with mis-matching ch's depending on method
(nrow(mismatch_ch) / nrow(ch_compare_df)) * 100  # percent of all fish with mis-matching ch's
  
# choose data set to proceed with
ch_df <- rk_ch    
    
# add metadata before saving
# note: this needs to be updated if we want to include covariates for supplementation fish and/or other adults that don't make it to dabom.
tag_summary <- ch_df %>%
  left_join(tag_meta %>%
              distinct(tag_code, .keep_all = TRUE) %>%
              # covariates to keep, for now
              select(spawn_yr,
                     tag_code,
                     srr,
                     spawn_site,
                     final_node,
                     mpg,
                     popid,
                     popname,
                     fl_mm = lgdf_lmm,
                     a_or_b,
                     fw_age,
                     sw_age,
                     total_age,
                     gen_sex,
                     gen_stock,
                     gen_stock_prob,
                     gen_parent_hatchery,
                     gen_by),
            by = c("spawn_yr", "tag_code"))
  # group_by(tag_code) %>%
  # filter(!(tag_code == "384.3B23AD5B2A" & row_number() > 1)) %>%
  # ungroup() %>%
  # mutate(across(c(SRR, LGDFLmm, GenSex, GenStock), ~ if_else(tag_code == "384.3B23AD5B2A", NA, .)))

# create model data
ch_mod_dat <- tag_summary %>%
  #select(-other) %>% # extra field from previous step
  filter(spawner_above == 1) %>%
  #filter(!(POP_NAME %in% c('GRA', 'RAPH', 'OXBO', 'DNFH'))) %>% #hatcheries and unknown spawning areas
  rowwise() %>%
  mutate(kelt_below = max(c_across(c(kelt_goa:rs_above))),
         kelt_rs = max(c_across(c(kelt_bon, rs_bon, rs_lgr, rs_above))),
         rs_all = max(c_across(c(rs_bon, rs_lgr, rs_above)))) %>%
  mutate(rs_lgr = max(c_across(c(rs_lgr, rs_above)))) %>%
  ungroup() %>%
  select(spawn_yr, tag_code, srr:gen_by, everything())

# a final summary of capture histories, by spawn year
ch_summary <- ch_mod_dat %>%
  group_by(spawn_yr) %>%
  summarise(across(.cols = c(release_lgr:rs_lgr, ), sum))

# save results
save(ch_mod_dat, file = paste0('./data/input/cjs_model_data_sy', max_sy, '.Rda'))

# clear environment
rm(list = ls())

### END SCRIPT