# Purpose: Load data from SnakeBasinStatus project and develop datasets for
# steelhead spawner, kelt and repeat spawner summaries.
# Author: Ryan N. Kinzer
# Date: 2021-10-08


# load pkgs----
library(tidyverse)
library(sf)
#library(lubridate)
#library(PITcleanr)
source('./R/identifyFishType.R')

# set script parameters----
reviewed <- TRUE
spp <- 'Steelhead'
yr_range <- 2010:2024
sy <- max(yr_range)
snake_proj <- '../SnakeRiverFishStatus/'
PITcleanr_folder <- paste0(snake_proj, 'output/PITcleanr/human_reviewed/')
Lifehistory_folder <- paste0(snake_proj, 'output/life_history/')
trap_path = paste0(snake_proj,
                   'data/LGTrappingDB/LGTrappingDB_2024-12-26.csv')

# load data----
load(file = paste0(snake_proj, 'data/configuration_files/site_config_LGR_20241226.rda')) ; rm(flowlines, parent_child, configuration, sr_site_pops)

clean_obs = list.files(path = PITcleanr_folder,
                          pattern = "Steelhead",
                          full.names = TRUE) %>%
  map_df(~ {
    # extract the spawn year from the file name
    spawn_yr <- str_extract(.x, "(?<=SY)\\d{4}")
    readxl::read_xlsx(path = .x, sheet = "Sheet1") %>%
      mutate(spawn_yr = as.numeric(spawn_yr))
  })

tag_meta <- list.files(path = Lifehistory_folder,
                       pattern = "Steelhead",
                       full.names = TRUE) %>%
  map_df(~ readxl::read_xlsx(
    path = .x,
    sheet = "tag_lh"
  )) %>%
  rename(spawn_yr = spawn_year)


trap_df = read_csv(trap_path) %>%
  mutate(GenStock = recode(GenStock, "LOWSALM" = "LOSALM"))

# get spawning location and subsequent kelt & repeat spawning observations
dat <- clean_obs %>%
  filter(!(life_stage == 'spawner' & user_keep_obs == FALSE)) # remove all false spawning locations.

# turn all obs prior to last spawner obs to spawner obs.
find_max_spawner <- function(df){
  tmp <- df %>%
    filter(life_stage == 'spawner') %>%
    group_by(tag_code) %>%
    slice(which.max(min_det)) %>%
    select(tag_code, max_spawner_det = min_det) %>%
    right_join(df) %>%
    mutate(tmp = case_when(
    min_det <= max_spawner_det  ~ 'spawner',
    TRUE ~ life_stage)
      )
  
  return(tmp)
    
}

# corrects observations of fish moving downstream looking similar to a kelt and then later move upstream following spawning patterns
dat <- dat %>%
  group_by(spawn_yr) %>%
  nest() %>%
  mutate(new = map(
    data,
    .f = ~find_max_spawner(.)
  )) %>%
  select(-data) %>%
  unnest(new) 

# Check QA/QC
# tmp <- new_obs %>%
#   filter(life_stage != tmp)

cth_df <- dat %>%
  mutate(life_stage = tmp) %>%
  select(-tmp) %>%
  mutate(site_code = str_remove(node, "_[D|U]")) %>%
  left_join(crb_sites_sf %>% select(site_code, rkm, rkm_total) %>% st_drop_geometry(), by = "site_code")


#extra to QA/QC
# find_post_kelt <- function(df){
#   tmp <- df %>%
#     filter(life_stage == 'kelt') %>%
#     group_by(tag_code) %>%
#     slice(which.min(min_det)) %>%
#     select(tag_code, min_kelt_det = min_det) %>%
#     right_join(df) %>%
#     filter(min_det >= min_kelt_det)
# 
#   return(tmp)
# }
# 
# post_kelt_obs <- new_obs %>%
#   group_by(spawn_yr) %>%
#   nest() %>%
#   mutate(filtered = map(
#     data,
#     .f = ~find_post_kelt(.)
#   )) %>%
#   select(-data) %>%
#   unnest(filtered)
# 
# tmp2 <- post_kelt_obs %>% group_by(spawn_yr, life_stage) %>% count()
# 
# 
# tags <- post_kelt_obs %>%
#   filter(life_stage == 'spawner') %>%
#   distinct(tag_code)
# 
# tmp3 <- inner_join(tmp, tags)
# 
# rm(tags, tmp, tmp2, tmp3)


# build capture-history matrix
ch_all <- cth_df %>%
  #filter(life_stage != "repeat spawner") %>%
  group_by(spawn_yr, tag_code) %>%
  mutate(lgr_max_det = max(max_det[site_code == "LGR"])) %>%
  group_by(spawn_yr) %>%
  nest() %>%
  mutate(ch = map(data,
                  .f = ~ steelheadCapHist(.))) %>%
  select(-data) %>%
  unnest(ch)
  
# Check QA/QC
tmp <- ch_all %>%
  group_by(spawn_yr) %>%
  summarise(across(-tag_code,sum))

# add covariate data----
tag_summary <- full_join(tag_meta, ch_all,
                         by = c('spawn_yr', 'tag_code'))


save(cth_df, tag_summary, file = paste0('./data/kelt_data_sy',sy,'.Rda'))



library(marked)

n = 1000
phi1 <- .3
phi2 <- .5
p1 <- .1
p2 <- .75

z1 = rep(1, n)
sum(z1)
z2 = rbinom(n, 1, phi1)
sum(z2)
z3 = z2 * rbinom(n, 1, phi2)
sum(z3)

y2 <- z2 * rbinom(n, 1, p1)
y3 <- z3 * rbinom(n, 1, p2)


df <- tibble(ch = paste0(z1, y2, y3))


ch_kelt <- ch_all %>%
  filter(spawner_above == 1) %>%
  mutate(ch = paste0(spawner_above, kelt_LGR, kelt_below)) %>%
  select(spawn_yr, ch) %>%
  ungroup()


# Convert data to marked format
m_data <- process.data(ch_kelt, model="CJS", groups = "spawn_yr")
m_data$begin.time <- 1
m_data$initial.age <- 0
# Add spawn year as a numeric covariate for Phi (survival)
design_data$Phi$spawn_yr <- factor(design_data$Phi$group)

# Add spawn year as a covariate for p (capture probability)
design_data$p$spawn_yr <- factor(design_data$p$group)

# Add time as a factor to allow time-varying effects
design_data$Phi$time <- factor(design_data$Phi$time)
design_data$p$time <- factor(design_data$p$time)

# Fit a Cormack-Jolly-Seber model (constant survival and capture probabilities)
cjs_model <- crm(m_data, model="CJS", hessian=TRUE, model.parameters=list(Phi=list(formula=~time + spawn_yr), p = list(formula=~time + spawn_yr)))

cjs_model$results

# Extract logit-scale estimates

Phi_hat <- exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$Phi)/(1 + exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$Phi))
p_hat <- exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$p)/(1 + exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$p))







# MA's summary

# kelt simple capture histories (new version); need to re-think what is actually a kelt obs at grs
kelt_ch_df = cth_df %>%
  filter(life_stage != "repeat spawner") %>%
  group_by(spawn_yr, tag_code) %>%
  mutate(lgr_max_det = max(max_det[site_code == "LGR"])) %>%
  summarise(
    spawner_above = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det & life_stage == 'spawner'), 1, 0),
    kelt_above = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det & life_stage == 'kelt'), 1, 0),
    #ups = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det), 1, 0), # RK commented
    kelt_grs = if_else(any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
    kelt_dwn = if_else(any(site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON") & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
    # when was the fish last observed moving upstream at lgr?
    lgr_max_det = unique(lgr_max_det),
    # if grs = 1, when was the fish last observed at grs as a kelt?
    grs_kelt_det = if (any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det)) {
      max(max_det[site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det], na.rm = TRUE)
    } else { NA },
    .groups = "drop"
  ) %>%
  # join some additional data from lgtrappingdb
  left_join(trap_df %>%
              filter(LGDLifeStage == "RF" & str_starts(SRR, "3")) %>%
              mutate(SpawnYear = as.numeric(str_remove(SpawnYear, "^SY"))) %>%
              select(spawn_yr = SpawnYear,
                     tag_code = LGDNumPIT,
                     collection_date = CollectionDate,
                     srr = SRR,
                     fl_mm = LGDFLmm,
                     gen_sex = GenSex,
                     gen_stock = GenStock) %>%
              group_by(spawn_yr, tag_code) %>%
              filter(collection_date == max(collection_date, na.rm = T)) %>%
              ungroup(),
            by = c("spawn_yr", "tag_code")) %>%
  # one tag code has errant data from lgtrappingdb; appears to be two separate fish sampled, but assigned to same tag_code
  group_by(tag_code) %>%
  filter(!(tag_code == "384.3B23AD5B2A" & row_number() > 1)) %>%
  ungroup() %>%
  mutate(across(c(srr, fl_mm, gen_sex, gen_stock), ~ if_else(tag_code == "384.3B23AD5B2A", NA, .))) %>%
  # remove hatchery fish
  filter(!str_ends(srr, "H")) %>%
  # join some additional data from lh_df
  left_join(tag_meta %>%
              distinct(tag_code, .keep_all = TRUE) %>%
              select(spawn_yr,
                     tag_code,
                     spawn_site,
                     mpg,
                     popid,
                     week_num),
            by = c("spawn_yr", "tag_code")) %>% # RK added distinct for a duplicate tag code: 3D9.1C2D0C3F01
  # calculate the julian week and month of when fish last left lgr
  mutate(julian_week = ceiling(yday(lgr_max_det) / 7)) %>%
  mutate(julian_month = month(lgr_max_det))


# check RK's vs. MA's

rk <- ch_all %>%
       filter(spawn_yr == 2010) %>%
       filter(spawner_above == 1) %>%
       pull(tag_code)

ma <- kelt_ch_df %>%
  filter(spawn_yr == 2010) %>%
  filter(spawner_above == 1) %>%
  pull(tag_code)
  
diff_codes <- rk[!(rk %in% ma)]

diff_obs <- cth_df %>%
  filter(spawn_yr == 2010) %>%
  filter(tag_code %in% diff_codes)



alpha <- 1 # priors for success and failures
beta <- 10

kelt_tbl = kelt_ch_df %>%
  #filter(pop_spwn == 1) %>%
  group_by(spawn_yr) %>%
  summarise(
    n_tags = n_distinct(tag_code),
    S = sum(spawner_above),
    #n_kelt_above = sum(kelt_above),
    t = sum(spawner_above == 1 & kelt_grs == 1),
    #n_kelt_grs = sum(kelt_grs == 1),
    n = sum(kelt_dwn == 1),
    x = sum(kelt_grs == 1 & kelt_dwn == 1),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  # left_join(valid_tag_df, by = "spawn_yr") %>%
  # select(spawn_yr,
  #        n_tags,
  #        everything()) %>%
  mutate(
    p = (x + alpha) / (n + alpha + beta),  # Bayesian smoothing
    #p_se = sqrt((p * (1 - p)) / n),
    k = t / p,
    #k_se = (t/p^2) * p_se,
    rate = pmin(k / S, 1),
    #rate_se = k_se/S,
    #lwr = pmax(rate - 1.96*rate_se, 0),
    #upr = pmin(rate + 1.96*rate_se, 1)
    )
# select(spawn_yr, n_tags, n, everything())

# Load necessary library
library(boot)

# Bootstrap function (requires 'data' and 'indices' arguments)
bootstrap_function <- function(data, indices) {
  # Extract values
  n <- data$n
  x <- data$x
  t <- data$t
  S <- data$S
  
  # Compute observed p
  p_obs <- x / n
  p_obs <- pmax(p_obs, 0.01)  # Avoid very small values
  
  # Resample x from Binomial(n, p_obs)
  x_boot <- rbinom(1, size = n, prob = p_obs)
  
  # Compute p, k_hat, and rate
  p_boot <- x_boot / n
  p_boot <- pmax(p_boot, 0.01)  # Avoid division by zero
  
  k_hat_boot <- t / p_boot
  rate_boot <- k_hat_boot / S
  rate_boot <- pmin(rate_boot, 1)  # Cap rate at 1.0
  
  return(rate_boot)
}

# Run bootstrap
set.seed(123) # For reproducibility

boot_results <- kelt_tbl %>%
  group_by(spawn_yr) %>%
  nest() %>%
  mutate(
    boot = map(data, ~ boot(data = ., statistic = bootstrap_function, R = 10000)),
    ci = map(boot, ~ tryCatch(boot.ci(., type = "perc", na.rm = TRUE), error = function(e) NULL)),
    rate_lwr = map_dbl(ci, ~ if (!is.null(.x)) .x$percent[4] else NA_real_),
    rate_upr = map_dbl(ci, ~ if (!is.null(.x)) .x$percent[5] else NA_real_)) %>%
  select(spawn_yr, rate_lwr, rate_upr)

# Merge results back into `kelt_tbl`
kelt_tbl <- kelt_tbl %>%
  left_join(boot_results, by = "spawn_yr")

 ggplot(data = kelt_tbl) +
   geom_ribbon(aes(x = spawn_yr, ymin = rate_lwr, ymax = rate_upr), alpha = .5) +
   geom_line(aes(x = spawn_yr, y = rate))
 