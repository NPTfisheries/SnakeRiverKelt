
#---------------
# MA's method for preparing data
# kelt capture histories
kelt_ch_df = dat %>%
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
              distinct(tag_code, .keep_all = TRUE) %>%
              filter(collection_date == max(collection_date, na.rm = T)) %>%
              ungroup(),
            by = c("spawn_yr", "tag_code")) %>%
  # one tag code has errant data from lgtrappingdb; appears to be two separate fish sampled, but assigned to same tag_code
  group_by(tag_code) %>%
  filter(!(tag_code == "384.3B23AD5B2A" & row_number() > 1)) %>%
  ungroup() %>%
  mutate(across(c(srr, fl_mm, gen_sex, gen_stock), ~ if_else(tag_code == "384.3B23AD5B2A", NA, .))) %>%
  # remove hatchery fish
  #filter(!str_ends(srr, "H")) %>%
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


#---------------
# compare RK & MA results
# sy = 2010
# rk = ch_all %>%
#   #filter(spawn_yr == sy) %>%
#   filter(spawner_above == 1) %>%
#   pull(tag_code)
# 
# ma = kelt_ch_df %>%
#   #filter(spawn_yr == sy) %>%
#   filter(spawn_abv == 1) %>%
#   pull(tag_code)
# 
# diff_codes <- rk[!(rk %in% ma)]
# 
# chk_h = trap_df %>%
#   filter(LGDNumPIT %in% diff_codes) 
# 
# diff_obs <- cth_df %>%
#   filter(spawn_yr == 2010) %>%
#   filter(tag_code %in% diff_codes)

# library(marked)
# 
# n = 1000
# phi1 <- .3
# phi2 <- .5
# p1 <- .1
# p2 <- .75
# 
# z1 = rep(1, n)
# sum(z1)
# z2 = rbinom(n, 1, phi1)
# sum(z2)
# z3 = z2 * rbinom(n, 1, phi2)
# sum(z3)
# 
# y2 <- z2 * rbinom(n, 1, p1)
# y3 <- z3 * rbinom(n, 1, p2)
# 
# df <- tibble(ch = paste0(z1, y2, y3))
# 
# ch_kelt <- ch_all %>%
#   filter(spawner_above == 1) %>%
#   mutate(ch = paste0(spawner_above, kelt_LGR, kelt_below)) %>%
#   select(spawn_yr, ch) %>%
#   ungroup()
# 
# # Convert data to marked format
# m_data <- process.data(ch_kelt, model="CJS", groups = "spawn_yr")
# m_data$begin.time <- 1
# m_data$initial.age <- 0
# # Add spawn year as a numeric covariate for Phi (survival)
# design_data$Phi$spawn_yr <- factor(design_data$Phi$group)
# 
# # Add spawn year as a covariate for p (capture probability)
# design_data$p$spawn_yr <- factor(design_data$p$group)
# 
# # Add time as a factor to allow time-varying effects
# design_data$Phi$time <- factor(design_data$Phi$time)
# design_data$p$time <- factor(design_data$p$time)
# 
# # Fit a Cormack-Jolly-Seber model (constant survival and capture probabilities)
# cjs_model <- crm(m_data, model="CJS", hessian=TRUE, model.parameters=list(Phi=list(formula=~time + spawn_yr), p = list(formula=~time + spawn_yr)))
# 
# cjs_model$results
# 
# # Extract logit-scale estimates
# 
# Phi_hat <- exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$Phi)/(1 + exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$Phi))
# p_hat <- exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$p)/(1 + exp(matrix(c(1,1,0,1), ncol = 2) %*% cjs_model$results$beta$p))

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
 