# Purpose: Estimate natural-origin kelt abundance back to GRS, GOA, and LMN.
#
# Author: Mike Ackerman & Ryan N. Kinzer
# Date Created: January 13, 2026
#   Date Last Modified: January 14, 2026

# load packages
library(tidyverse)
library(marked)
library(readxl)

# set arguments
max_sy = 2024

# load data from 01_prep_data.R
load(paste0("./data/input/compiled_lgr_kelt_data_sy", max_sy, ".rda"))

mod_df = ch_bio_df %>%
  # focus on fish with spawning observations and a final spawning population (note: this excludes nearly all hatchery fish bc they didn't make it to dabom)
  #filter(spawner_above == 1 & !is.na(popid)) %>%
  # filter out a few remaining hatchery fish (i.e., focus on known 32W)
  filter(srr == "32W") %>%
  # recode some spawners to populations
  mutate(popid = recode(popid, "CRLMA-s/CRSFC-s" = "CRSFC-s")) %>%
  # remove remaining spawners with indeterminant spawning population (SFG, USE, USI)
  filter(!grepl("/", popid)) %>%
  # recode some spawners to populations
  # mutate(popid = case_when(
  #   popid == "CRLMA-s/CRSFC-s"         ~ "CRSFC-s",
  #   popid == "SFMAI-s/SFSEC-s"         ~ "und_SFMAI-s",
  #   popid == "SRPAH-s/SREFS-s/SRUMA-s" ~ "und_SRPAH-s",
  #   TRUE ~ popid
  # )) %>%
  # trim to fish with known sex
  filter(gen_sex %in% c("F", "M")) %>%
  # finally, remove kelts that were removed by kelt reconditioning program
  filter(reconditioned_kelt == "not_reconditioned") %>%
  rowwise() %>%
  mutate(kelt_blw = max(c_across(c(kelt_ihr:rs_above))),                   # was kelt observed anywhere downstream, including on return spawn
         kelt_rs  = max(c_across(c(kelt_bon, rs_bon, rs_lgr, rs_above))),  # was kelt observed making it down to bon and/or as return spawner
         rs_all   = max(c_across(c(rs_bon, rs_lgr, rs_above)))) %>%        # was kelt observed as a repeat spawner to bon or upstream
  mutate(rs_lgr   = max(c_across(c(rs_lgr, rs_above)))) %>%                # was kelt observed as a repeat spawner to lgr or upstream
  ungroup() %>%
  select(spawn_yr, tag_code, srr:reconditioned_kelt, everything())

# review capture histories, by spawn year
mod_df %>%
  group_by(spawn_yr) %>%
  summarise(release_lgr = sum(release_lgr),
            spawner_above = sum(spawner_above),
            kelt_grs = sum(kelt_grs),
            kelt_goa = sum(kelt_goa),
            kelt_lma = sum(kelt_lma),
            kelt_blw = sum(kelt_blw))

#-----------------------------------------------------------------------
# Estimate total natural-origin kelt abundance back to GRS, GOA, and LMN

#---------------
# set up and run cjs models

mod_dat = data.frame(
  ch = paste0(mod_df$release_lgr, mod_df$kelt_grs, mod_df$kelt_goa, mod_df$kelt_lma, mod_df$kelt_blw),
  freq = rep(1, dim(mod_df)[1]),
  #mpg = as.factor(mod_df$mpg),
  #pop = as.factor(mod_df$popid),
  year = as.factor(mod_df$spawn_yr),
  sex = as.factor(mod_df$gen_sex)
)

# create function to run all models
fit.kelt.cjs.models <- function(){
  
  # Apparent survival (Phi) formulas:
  
  # Baseline models
  #Phi.dot          <- list(formula = ~ 1)   # constant survival: all kelts have same probability of surviving to LGR, regardless of year, MPG, population, or sex
  Phi.time         <- list(formula = ~ time) # survival to LGR differs across kelts, but is constant across years, MPG, population, and sex
  
  # Additive Models
  Phi.sex          <- list(formula = ~ sex + time)  # survival to LGR varies by sex (consistent across mpg/pop, and year)
  #Phi.mpg          <- list(formula = ~ mpg + time)  # survival to LGR varies by mpg (consistent across sex and year)        
  #Phi.pop          <- list(formula = ~ pop + time)  # survival to LGR varies by pop (consistent across sex and year) 
  Phi.year         <- list(formula = ~ year + time) # survival to LGR varies by year (consistent across mpg/pop, and sex)
  
  # Additive Effects (Multiple Covariates)
 # Phi.pop.sex      <- list(formula = ~ pop + sex + time)  # survival varies by pop and sex (consistent across year)
  #Phi.mpg.sex      <- list(formula = ~ mpg + sex + time)  # survival varies by mpg and sex (consistent across year)
  #Phi.year.pop     <- list(formula = ~ year + pop + time) # survival varies by year and pop (consistent across sex)
  #Phi.year.mpg     <- list(formula = ~ year + mpg + time) # survival varies by year and mpg (consistent across sex)
  Phi.year.sex     <- list(formula = ~ year + sex + time) # survival varies by year and sex (consistent across mpg/pop)
  
  #Phi.year.mpg.sex <- list(formula = ~ year + mpg + sex + time) # survival varies by year, mpg, and sex
  #Phi.year.pop.sex <- list(formula = ~ year + pop + sex + time) # survival varies by year, pop, and sex
  #Phi.year.pop.sex.int <- list(formula = ~ year + pop * sex + time)
  #Phi.mpg.pop.sex  <- list(formula = ~ mpg + pop + sex + time)  # survival varies by mpg, pop, and sex (consistent across year)
  
  # Interactions
  #Phi.mpg.time     <- list(formula = ~ mpg * time)        # survival differences among MPGs change over time
  #Phi.pop.time     <- list(formula = ~ pop * time)        # survival differences among populations change over time
  Phi.year.time    <- list(formula = ~ year * time)       # survival differences across years change over time
  Phi.sex.time     <- list(formula = ~ sex * time)        # survival differences between sexes change over time
  
  # Fully Interactive
  #Phi.full         <- list(formula = ~ year * pop * sex * time) # (mpg removed) survival to LGR is affected by year, pop, and sex, with all factors interacting

  # Detection probability (p) formulas
  #p.pop.time      <- list(formula = ~ pop * time) 
  p.year.time     <- list(formula = ~ year * time)
  p.year.sex.time <- list(formula = ~ year * time + sex)
  #p.sex.time      <- list(formula = ~ sex * time)
  #p.pop           <- list(formula = ~ pop + time)
  #p.year.sex      <- list(formula = ~ year + sex + time)
  #p.year          <- list(formula = ~ year + time)
  #p.sex           <- list(formula = ~ sex + time)
  #p.time          <- list(formula = ~ time)        # one discrete estimate of p per capture event
  #p.dot           <- list(formula = ~ 1)           # constant detection
  
  # Construct all combinations and put into one model table
  cml     <- create.model.list(c("Phi","p")) # makes all possibile combinations of those parameter formulas
  results <- crm.wrapper(cml, 
                         data = kelt.proc, 
                         ddl = kelt.ddl,
                         external = FALSE, 
                         accumulate = FALSE, 
                         hessian = TRUE)
  return(results)
}

# Phi.all <- list(formula=~pop*year*time)
# p.all <- list(formula=~pop*year*time)
# 
# mod.all <- crm(kelt.proc, kelt.ddl,
#                model.parameters = list(Phi = Phi.all,
#                                        p = p.all),
#                accumulate = FALSE)

kelt.proc <- process.data(mod_dat, begin.time = 1, model = 'CJS', groups = c("year", "sex"))
kelt.ddl  <- make.design.data(kelt.proc)

kelt.cjs.models <- fit.kelt.cjs.models()

kelt.cjs.models
best_mod = kelt.cjs.models[[9]] # Phi(~year + sex + time)p(~year * time + sex)

best_res <- best_mod$results$reals %>%
  bind_rows(.id = 'param') %>%
  mutate(
    reach = case_when(
      param == "Phi" & occ == 1 ~ "LGR to GRS",
      param == "Phi" & occ == 2 ~ "GRS to GOA",
      param == "Phi" & occ == 3 ~ "GOA to LMA",
      param == "Phi" & occ == 4 ~ "LMA to Down",
      param == "p"   & occ == 2 ~ "GRS",
      param == "p"   & occ == 3 ~ "GOA",
      param == "p"   & occ == 4 ~ "LMA",
    )
  ) %>%
  filter(!occ == 5)

lgr_N_df = read_xlsx(path = "./data/input/Lower Granite Abundance Estimates.xlsx",
                   sheet = "Wild_Sex") %>%
  filter(SpeciesName == "Steelhead",
         Sex %in% c("Female", "Male"),
         SpawnYear %in% 2010:2024) %>%
  mutate(
    Sex = case_when(
      Sex == "Female" ~ "F",
      Sex == "Male"   ~ "M",
      TRUE            ~ Sex
    )
  )

reach_levels = c("LGR to GRS", "GRS to GOA", "GOA to LMA", "LMA to Down")

kelt_esc_df = best_res %>%
  filter(param == "Phi") %>%
  transmute(
    year,
    sex,
    reach,
    phi = estimate,
    location = case_when(
      reach == "LGR to GRS"  ~ "GRS",
      reach == "GRS to GOA"  ~ "GOA",
      reach == "GOA to LMA"  ~ "LMA",
      reach == "LMA to Down" ~ "Down",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(reach = factor(reach, levels = reach_levels)) %>%
  arrange(year, sex, reach) %>%
  group_by(year, sex) %>%
  mutate(cum_phi = cumprod(phi)) %>%
  ungroup() %>%
  left_join(lgr_N_df %>%
              transmute(year = as.factor(SpawnYear), sex = Sex, esc_lgr = Estimate),
            by = c("year", "sex")) %>%
  mutate(escapement = esc_lgr * cum_phi) %>%
  #filter(reach == "LGR to GRS") %>% select(-cum_phi, -location)
  select(year, sex, location, escapement) %>%
  bind_rows(lgr_N_df %>%
              transmute(year = as.factor(SpawnYear), sex = Sex, location = "LGR", escapement = Estimate)) %>%
  mutate(location = factor(location, levels = c("LGR", "GRS", "GOA", "LMA", "Down"))) %>%
  arrange(year, sex, location) %>%
  filter(sex == "F" & location != "Down") %>%
  mutate(escapement = round(escapement)) %>%
  pivot_wider(names_from = location,
              values_from = escapement)
 
# write results
write_csv(kelt_esc_df,
          file = "./output/kelt_abund_2_sr_dams.csv") 

#-------------------
# simple Monte-Carlo
# n_sims = 5000
# 
# # assume abundance is lognormal (get lognormal params from median & CI)
# ln_pars = function(med, lcl, ucl) {
#   sdlog   = (log(ucl) - log(lcl)) / (2 * 1.96)
#   meanlog = log(med)
#   list(meanlog = meanlog, sdlog = sdlog)
# }
# 
# lgr_N_draws = lgr_N_df %>%
#   rename(year = SpawnYear,
#          sex = Sex,
#          med = Estimate,
#          lcl = L_CI,
#          ucl = U_CI) %>%
#   mutate(pars = pmap(list(med, lcl, ucl), ln_pars)) %>%
#   mutate(draws = map(pars, ~ rlnorm(n_sims, .$meanlog, .$sdlog))) %>%
#   select(year, sex, draws)
# 
# # assume survival is logit-normal
# logit_pars = function(mu, lcl, ucl) {
#   logit <- function(p) log(p / (1 - p))
#   sd <- (logit(ucl) - logit(lcl)) / (2 * 1.96)
#   mean <- logit(mu)
#   list(mean = mean, sd = sd)
# }
# 
# # add epsilon to squish into (0,1)
# eps <- 1e-6
# 
# phi_draws = best_res %>%
#   filter(param == "Phi") %>%
#   transmute(
#     year, sex, reach,
#     mu  = pmin(pmax(estimate, eps), 1 - eps),
#     lcl = pmin(pmax(lcl,     eps), 1 - eps),
#     ucl = pmin(pmax(ucl,     eps), 1 - eps)
#   ) %>%
#   mutate(
#     pars  = pmap(list(mu, lcl, ucl), logit_pars),
#     draws = map(pars, ~ plogis(rnorm(n_sims, mean = .x$mean, sd = .x$sd)))
#   ) %>%
#   select(year, sex, reach, draws)
# 
# escape_draws = lgr_draws 

# save important results
#save(mod_df, mod_dat, kelt.cjs.models, best_mod, best_res, file = paste0('./data/output/mod_dat_sy', max_sy, '_', Sys.Date(), '.rda'))

### END SCRIPT
