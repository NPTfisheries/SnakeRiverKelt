# Purpose: Load data from SnakeRiverFishStatus project and develop inputs for
#   steelhead spawner, kelt, and repeat spawner summaries.
#
# Author: Ryan N. Kinzer
# Date Created: 2021-10-08
#   Date Last Modified: 2026-01-13
#   Modified By: Mike Ackerman

# load packages
library(tidyverse)
library(sf)
library(readxl)
library(janitor)

# source some helper functions
source('./R/find_max_spawner.R')
source("./R/steelheadCapHist.R")

# set some parameters
yr_range = 2010:2024
max_sy = max(yr_range)

#---------------
# load data

# set some file paths to datasets
snake_proj = "../SnakeRiverFishStatus/"
PITcleanr_folder = paste0(snake_proj, "output/PITcleanr/human_reviewed/")
lh_folder = paste0(snake_proj, "output/life_history/")
trap_path = paste0(snake_proj, "data/LGTrappingDB/LGTrappingDB_2026-01-06.csv")

# load configuration files (just need crb_sites_sf object)
load(file = paste0(snake_proj, "data/configuration_files/site_config_LGR_20250416.rda")) ; rm(flowlines, parent_child, configuration, sr_site_pops)

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

# load life history data (do I need this?)
lh_df = list.files(path = lh_folder,
                   pattern = "Steelhead",
                   full.names = T) %>%
  map_df(~ readxl::read_xlsx(
    path = .x,
    sheet = "tag_lh"
  )) %>%
  rename(spawn_yr = spawn_year) %>%
  janitor::clean_names()

# load lgtrappingdb
trap_df = read_csv(trap_path) %>%
  janitor::clean_names() %>%
  # some LOSALM & LOCLWR genstocks are errantly coded as LOWSALM & LOWCLWR
  mutate(gen_stock = recode(gen_stock, "LOWSALM" = "LOSALM")) %>%
  mutate(gen_stock = recode(gen_stock, "LOWCLWR" = "LOCLWR")) %>%
  # keep just adults (i.e., returning fish)
  filter(lgd_life_stage == "RF",
         !spawn_year == "None") %>%
  # clean up some columns
  rename(lgd_fl_mm = lgdf_lmm,
         spawn_yr = spawn_year,
         tag_code = lgd_num_pit) %>%
  mutate(spawn_yr = as.numeric(str_remove(spawn_yr, "SY")))

# load kelts collected for reconditioning program from critfc
critfc_kelt_df = read_xlsx(path = "data/input/SY2016-2021 collected fish with tag_lej.xlsx",
                           sheet = "16-21 SY existing tag list") %>%
  janitor::clean_names() %>%
  # clean up collection_year column
  mutate(collection_year = as.numeric(recode(collection_year, "2019 - marked SY2018" = "2019")),
         reconditioned_kelt = "existing_pit_tag") %>%
  # join info to parse kelts newly tagged at lgr vs. those with previous tags
  left_join(
    read_xlsx(path = "data/input/SY2016-2021 collected fish with tag_lej.xlsx",
              sheet = "16-21 LGR tagged and collected") %>%
      janitor::clean_names() %>%
      select(collection_year, pit_tag_code) %>%
      mutate(new_lgr_tag = TRUE),
    by = c("collection_year", "pit_tag_code")
  ) %>%
  mutate(reconditioned_kelt = if_else(!is.na(new_lgr_tag), "new_lgr_pit_tag", reconditioned_kelt)) %>%
  select(-new_lgr_tag,
         tag_code = pit_tag_code)

# quick qc; there's some minor mismatches from critfc's summary
tabyl(critfc_kelt_df, reconditioned_kelt, collection_year) %>%
  adorn_totals()

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
fix_errant_kelt_calls = TRUE

# apply find_max_spawner only if fix_errant_kelt_calls is TRUE
if (fix_errant_kelt_calls == TRUE) {
  cth_df = cth_df %>%
    group_by(spawn_yr) %>%
    nest() %>%
    mutate(new = map(data, find_max_spawner)) %>%
    select(-data) %>%
    unnest(new) 
}

# RK method
rk_ch = cth_df %>%
  ungroup() %>%
  group_by(spawn_yr) %>%
  nest() %>%
  mutate(ch = map(data,
                  .f = ~steelheadCapHist(.))) %>%
  select(-data) %>%
  unnest(ch)

# MA method
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
  mutate(user = "ma") %>%
  bind_rows(rk_ch %>%
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
ch_df = rk_ch

# add metadata
ch_bio_df = ch_df %>%
  # from lower granite trapping database
  left_join(trap_df %>%
              filter(!is.na(tag_code)) %>%
              distinct(spawn_yr, tag_code, .keep_all = TRUE) %>%
              select(spawn_yr,
                     tag_code,
                     collection_date,
                     srr,
                     lgd_marks_all,
                     lgd_mark_ad,
                     lgd_tags_all,
                     lgd_fl_mm,
                     bio_scale_final_age,
                     ptagis_last_event_site,
                     ptagis_last_event_date,
                     ptagis_event_last_spawn_site,
                     gen_sex,
                     gen_stock,
                     #gen_stock_prob,
                     gen_parent_hatchery,
                     gen_by,
                     gen_pbt_by_hat,
                     gen_pbt_r_group),
            by = c("spawn_yr", "tag_code")) %>%
  # additional fields from SnakeRiverFishStatus life history folders (this will only include fish that made it to dabom; e.g., exclude supplementation fish)
  left_join(lh_df %>%
              distinct(spawn_yr, tag_code, .keep_all = TRUE) %>%
              select(spawn_yr,
                     tag_code,
                     spawn_site,
                     final_node,
                     mpg,
                     popid,
                     popname,
                     fw_age,
                     sw_age,
                     total_age),
            by = c("spawn_yr", "tag_code")) %>%
  # identify fish within the kelt reconditioning program dataset
  left_join(critfc_kelt_df %>%
              select(collection_year, tag_code, reconditioned_kelt),
            by = c("spawn_yr" = "collection_year", "tag_code")) %>%
  mutate(reconditioned_kelt = ifelse(is.na(reconditioned_kelt), "not_reconditioned", reconditioned_kelt))

#---------------
# some qa/qc

# are there any duplicate tag codes within a spawn year?
any(duplicated(ch_bio_df[c("spawn_yr", "tag_code")]))

# check capture histories of reconditioned kelts
ch_bio_df %>%
  filter(!reconditioned_kelt == "not_reconditioned") %>%
  select(spawn_yr,
         tag_code,
         release_lgr:rs_above) %>%
  group_by(spawn_yr) %>%
  summarise(across(-tag_code, sum)) %>%
  ungroup()

# check sexes of reconditioned vs. all kelts
ch_bio_df %>%
  tabyl(reconditioned_kelt, gen_sex)

# check srrs of reconditioned vs. all kelts
ch_bio_df %>%
  tabyl(reconditioned_kelt, srr)

# data frame of 32H reconditioned kelts
rc_kelt_32h = ch_bio_df %>%
  filter(!reconditioned_kelt == "not_reconditioned",
         srr == "32H")

#---------------
# a final summary of capture histories, by spawn year
ch_summary = ch_bio_df %>%
  group_by(spawn_yr) %>%
  summarise(across(.cols = c(release_lgr:rs_above, ), sum))

# save results
save(ch_bio_df, file = paste0("./data/input/compiled_lgr_kelt_data_sy", max_sy, ".rda"))

# clear environment
rm(list = ls())

### END SCRIPT