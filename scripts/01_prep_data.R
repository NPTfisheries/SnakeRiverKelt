# Purpose: Load data from SnakeRiverFishStatus project and develop inputs for
#   steelhead spawner, kelt, and repeat spawner summaries.
#
# Author: Ryan N. Kinzer
# Date Created: 2021-10-08
#   Date Modified: 2025-03-13
#   Modified By: Mike Ackerman

# load packages and functions
library(tidyverse)
library(sf)
#library(janitor)

source('./R/find_max_spawner.R')
source("./R/steelheadCapHist.R")

# set some parameters
#reviewed = TRUE
spp = "Steelhead"
yr_range = 2010:2024
sy = max(yr_range)
snake_proj = "../SnakeRiverFishStatus/"
PITcleanr_folder = paste0(snake_proj, "output/PITcleanr/human_reviewed/")
lh_folder = paste0(snake_proj, "output/life_history/")
trap_path = paste0(snake_proj, "data/LGTrappingDB/LGTrappingDB_2024-12-26.csv")

# load configuration files
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
  rename(spawn_yr = spawn_year)

# load lgtrappingdb
trap_df = read_csv(trap_path) %>%
  mutate(GenStock = recode(GenStock, "LOWSALM" = "LOSALM"))

#---------------
# Generate a common dataset

dat = pitcleanr_df %>%
  # remove spawner observations considered FALSE for dabom
  filter(!(life_stage == "spawner" & user_keep_obs == FALSE)) %>%
  mutate(site_code = str_remove(node, "_[D|U]")) %>%
  # correct obs of fish moving downstream looking similar to a kelt, then later moving upstream following spawning patters
  group_by(spawn_yr) %>%
  nest() %>%
  mutate(new = map(
    data,
    .f = ~find_max_spawner(.)
  )) %>%
  select(-data) %>%
  unnest(new) %>%
  group_by(spawn_yr, tag_code) %>%
  mutate(lgr_max_det = max(max_det[site_code == "LGR" & life_stage == "spawner"])) %>%
  ungroup() %>%
  select(spawn_yr,
         id,
         tag_code,
         life_stage,
         site_code,
         node,
         direction,
         slot,
         min_det,
         max_det,
         max_spawner_det,
         lgr_max_det,
         path) %>%
  # join rkm for each site
  left_join(crb_sites_sf %>%
              select(site_code,
                     rkm,
                     rkm_total) %>%
              st_drop_geometry(),
            by = "site_code")

#RK direction

rk_all <- dat %>%
  ungroup() %>%
  group_by(spawn_yr) %>%
  nest() %>%
  mutate(ch = map(data,
                  .f = ~steelheadCapHist(.))) %>%
  select(-data) %>%
  unnest(ch)

#MK direction

ma_all = dat %>%
# remove repeat spawner observations, for now
#filter(life_stage != "repeat spawner") %>%
group_by(spawn_yr, tag_code) %>%
summarise(
  release_lgr = if_else(any(max_det == lgr_max_det & life_stage == "spawner" & node == "LGR"), 1, 0),
  spawner_above = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det & life_stage == "spawner"), 1, 0),
  spawner_below = if_else(any(life_stage == "spawner" & grepl("GRS", path) & min_det > lgr_max_det & !site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON")), 1, 0),
  kelt_above = if_else(any(life_stage == "kelt"    & !grepl("GRS", path) & node != "LGR"), 1, 0),
  kelt_grs = if_else(any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  kelt_goa = if_else(any(site_code == "GOA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  kelt_lma = if_else(any(site_code == "LMA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  kelt_ihr = if_else(any(site_code == "IHR" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  kelt_mcn = if_else(any(site_code == "MCN" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  kelt_jda = if_else(any(site_code == "JDA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  kelt_tda = if_else(any(site_code == "TDA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  kelt_bon = if_else(any(site_code == "BON" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  #kelt_dwn = if_else(any(site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON") & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
  rs_bon   = if_else(any(life_stage == "repeat spawner" & node == "BON"), 1, 0),
  rs_lgr   = if_else(any(life_stage == "repeat spawner" & node == "LGR"), 1, 0),
  rs_above = if_else(any(life_stage == "repeat spawner" & !grepl("GRS", path) & node != "LGR"), 1, 0),
  # when was the fish last observed moving upstream at lgr?
  lgr_max_det = unique(lgr_max_det),
  # if grs = 1, when was the fish last observed at grs as a kelt?
  grs_kelt_det = if (any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det)) {
    max(max_det[site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det], na.rm = TRUE)
  } else { NA },
  .groups = "drop"
)
    
  # QA/QC: summed cap hists by spawn year
ma_ch_summary = ma_all %>%
  select(spawn_yr, 
         tag_code,
         release_lgr:rs_above) %>%
  group_by(spawn_yr) %>%
  summarise(across(-tag_code, sum)) %>%
  ungroup()

# QA/QC: summed cap hists by spawn year
rk_ch_summary = rk_all %>%
  group_by(spawn_yr) %>%
  summarise(across(-tag_code, sum)) %>%
  ungroup() %>%
  select(-other)

# join ma & rk capture histories
ch_compare_df = ma_all %>%
  select(spawn_yr:rs_above) %>%
  mutate(user = "ma") %>%
  bind_rows(rk_all %>%
              select(-other) %>%
              mutate(user = "rk")) %>%
  select(spawn_yr, tag_code, user, everything()) %>%
  arrange(spawn_yr, tag_code, user)

# return records where capture histories don't match (takes awhile)
mismatches = ch_compare_df %>%
  group_by(spawn_yr, tag_code) %>%
  filter(any(across(release_lgr:rs_above, ~ . != first(.)))) %>%
  ungroup()
  
nrow(mismatches) / 2
  (nrow(mismatches) / nrow(ch_compare_df)) * 100
  
# Choose data set to proceed with
ch_all <- rk_all    
    
# Add metadata before saving
tag_summary <- ch_all %>%
  left_join(tag_meta %>%
              distinct(tag_code, .keep_all = TRUE),
            by = c("spawn_yr", "tag_code")) %>%
  group_by(tag_code) %>%
  filter(!(tag_code == "384.3B23AD5B2A" & row_number() > 1)) %>%
  ungroup() %>%
  mutate(across(c(SRR, LGDFLmm, GenSex, GenStock), ~ if_else(tag_code == "384.3B23AD5B2A", NA, .)))

# create model data----
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
  select(spawn_yr, tag_code, species:a_or_b, everything())

tmp <- ch_mod_dat %>%
  group_by(spawn_yr) %>%
  summarise(across(.cols = c(release_lgr:rs_lgr, ), sum))

save(ch_mod_dat, file = paste0('./data/input/cjs_model_data_sy',sy,'.Rda'))

rm(list = ls())

