# Purpose: Final prep of results from 01_prep_data for CJS modes
#   and run various CJS models
#
# Author: Ryan N. Kinzer
# Date Created: Unknown
#   Date Last Modified: 2025-04-02
#   Modified By: Mike Ackerman

# load packages
library(tidyverse)
library(marked)

# set arguments
max_sy = 2024

# load data from 01_prep_data.R
load(paste0("./data/input/compiled_lgr_kelt_data_sy", max_sy, ".rda"))

# final prep and filtering for cjs models
mod_df = ch_bio_df %>%
  # focus on fish with spawning observations and a final spawning population (note: this excludes nearly all hatchery fish bc they didn't make it to dabom)
  filter(spawner_above == 1 & !is.na(popid)) %>%
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
  mutate(kelt_blw = max(c_across(c(kelt_goa:rs_above))),                 # was kelt observed anywhere downstream, including on return spawn
         kelt_rs = max(c_across(c(kelt_bon, rs_bon, rs_lgr, rs_above))), # was kelt observed making it down to bon and/or as return spawner
         rs_all = max(c_across(c(rs_bon, rs_lgr, rs_above)))) %>%        # was kelt observed as a repeat spawner to bon or upstream
  mutate(rs_lgr = max(c_across(c(rs_lgr, rs_above)))) %>%                # was kelt observed as a repeat spawner to lgr or upstream
  ungroup() %>%
  select(spawn_yr, tag_code, srr:reconditioned_kelt, everything())

# review capture histories, by spawn year
mod_df %>%
  group_by(spawn_yr) %>%
  summarise(spawner_above = sum(spawner_above),
            kelt_grs = sum(kelt_grs),
            kelt_blw = sum(kelt_blw))

#---------------
# set up and run cjs models

mod_dat = data.frame(
  ch = paste0(mod_df$spawner_above, mod_df$kelt_grs, mod_df$kelt_blw),
  freq = rep(1, dim(mod_df)[1]),
  mpg = as.factor(mod_df$mpg),
  pop = as.factor(mod_df$popid),
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
  Phi.mpg          <- list(formula = ~ mpg + time)  # survival to LGR varies by mpg (consistent across sex and year)        
  Phi.pop          <- list(formula = ~ pop + time)  # survival to LGR varies by pop (consistent across sex and year) 
  Phi.year         <- list(formula = ~ year + time) # survival to LGR varies by year (consistent across mpg/pop, and sex)
  
  # Additive Effects (Multiple Covariates)
  Phi.pop.sex      <- list(formula = ~ pop + sex + time)  # survival varies by pop and sex (consistent across year)
  Phi.mpg.sex      <- list(formula = ~ mpg + sex + time)  # survival varies by mpg and sex (consistent across year)
  Phi.year.pop     <- list(formula = ~ year + pop + time) # survival varies by year and pop (consistent across sex)
  Phi.year.mpg     <- list(formula = ~ year + mpg + time) # survival varies by year and mpg (consistent across sex)
  Phi.year.sex     <- list(formula = ~ year + sex + time) # survival varies by year and sex (consistent across mpg/pop)
  
  Phi.year.mpg.sex <- list(formula = ~ year + mpg + sex + time) # survival varies by year, mpg, and sex
  Phi.year.pop.sex <- list(formula = ~ year + pop + sex + time) # survival varies by year, pop, and sex
  #Phi.mpg.pop.sex  <- list(formula = ~ mpg + pop + sex + time)  # survival varies by mpg, pop, and sex (consistent across year)
  
  # Interactions
  Phi.mpg.time     <- list(formula = ~ mpg * time)        # survival differences among MPGs change over time
  Phi.pop.time     <- list(formula = ~ pop * time)        # survival differences among populations change over time
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

kelt.proc <- process.data(mod_dat, begin.time = 1, model = 'CJS', groups = c("mpg", "pop", "year", "sex"))
kelt.ddl  <- make.design.data(kelt.proc)

kelt.cjs.models <- fit.kelt.cjs.models()

kelt.cjs.models
#best_mod <- kelt.cjs.models[[24]] #$Phi.pop.p.pop.time
best_mod = kelt.cjs.models[[28]] # Phi.year.pop.sex

best_res <- best_mod$results$reals %>%
  bind_rows(.id = 'param')

#---------------
# set up and run cjs models

# marginize over the year covariate to calcualte an average survival for each pop
#source('./R/theme_rk.R')

# apparent survival, by year and sex
fig_s = best_res %>%
  filter(time == 1) %>%
  filter(param == "Phi") %>%
  ggplot(aes(x = fct_reorder(pop, estimate), y = estimate, color = sex)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl),
                width = 0.5,
                position = position_dodge(width = 0.5)) +
  geom_point(             position = position_dodge(width = 0.5)) +
  coord_flip() +
  facet_wrap(~year, ncol = 5) +
  labs(title = "Snake River Basin Steelhead Kelt",
       subtitle = "Survival of adult steelhead from Snake River basin spawning areas to the kelt life-stage at Lower Granite Dam.",
       x = "",
       y = "P(Survival | Observed Spawner)",
       color = "Sex") +
  scale_color_manual(values = c("F" = "#E9967A", "M" = "#008080")) + 
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),        
        axis.text = element_text(size = 8),                           
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0),
        legend.position = "right")
fig_s  

# fig_s <- best_res %>%
#   filter(time == 1) %>%
#   filter(param == 'Phi') %>%
#   #filter(estimate < .99) %>%
#   ggplot(aes(x = fct_reorder(pop,estimate), y = estimate)) +
#   geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
#   geom_point() +
#   coord_flip() +
#   facet_wrap(~year, ncol= 5) +
#   labs(title = 'Snake River Basin Steelhead Kelt',
#        subtitle = 'Survival of adult steelhead from Snake River basin spawning areas to the kelt life-stage at Lower Granite Dam.',
#        x = '',
#        y = 'P(Survival | Observed Spawner)') +
#   theme_classic() +
#   theme(axis.title = element_text(size = 12, face = "bold"),          # Larger axis labels
#         axis.text = element_text(size = 8),                           # Readable tick labels
#         legend.title = element_text(size = 12, face = "bold"),
#         legend.text = element_text(size = 10),
#         plot.title = element_text(size = 12, face = "bold", hjust = 0)) # Centered title

ggsave('./figures/kelt_lgd_survival_sex.png', fig_s,
       width = 11, height = 8.5,
       dpi = 600, units = 'in', device = "png")

# annual kelt detection probability at lgr
fig_d = best_res %>%
  # time 2 is lgr
  filter(time == 2) %>%
  filter(param == 'p') %>%
  ggplot(aes(x = fct_rev(year), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.4) +
  geom_point() +
  coord_flip() +
  labs(title = 'Snake River Basin Steelhead Kelt',
       subtitle = 'Detection probability of steelhead kelt at Lower Granite Dam after being observed in Snake River basin spawning areas.',
       x = 'Spawn Year',
       y = 'P(Detection | Kelt at LGD)') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0))
fig_d

ggsave('./figures/kelt_lgd_detection.png', fig_d,
       width = 11, height = 8.5,
       dpi = 600, units = 'in', device = "png")

# pop average - uses an inverse weighted mean to account for precision around each individual estimate
phi_lgd = best_res %>%
  filter(time == 1) %>%
  filter(param == 'Phi') %>%
  group_by(pop, sex) %>%
  summarize(
    avg_Phi = sum(estimate / (se^2), na.rm = TRUE) / sum(1 / (se^2), na.rm = TRUE),
    se_Phi = sqrt(1 / sum(1 / (se^2), na.rm = TRUE)),  # Weighted SE
    lower_CI = pmax(avg_Phi - 1.96 * se_Phi, 0),       # 95% CI lower bound
    upper_CI = pmin(avg_Phi + 1.96 * se_Phi, 1),       # 95% CI upper bound
    .groups = "drop"
  )

# average population survival, by sex
fig_avg_s = phi_lgd %>%
  ggplot(aes(x = fct_reorder(pop, avg_Phi), y = avg_Phi, color = sex)) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.4) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = .25)) +
  geom_point() +
  coord_flip() +
  labs(title = 'Snake River Basin Steelhead Kelt',
       subtitle = 'Average survival of adult steelhead from Snake River population spawning areas to the kelt life-stage at Lower Granite Dam.',
       x = '',
       y = 'Average Survival',
       color = "Sex") +
  scale_color_manual(values = c("F" = "#E9967A", "M" = "#008080")) +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),  # Larger axis labels
        axis.text = element_text(size = 12),  # Readable tick labels
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0))
fig_avg_s
  
# fig_avg_s <- phi_lgd %>%
#   ggplot(aes(x = fct_reorder(pop,avg_Phi), y = avg_Phi)) +
#   geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI)) +
#   scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = .25)) +
#   geom_point() +
#   coord_flip() +
#   labs(title = 'Snake River Basin Steelhead Kelt',
#        subtitle = 'Average survival of adult steelhead from Snake River population spawning areas to the kelt life-stage at Lower Granite Dam.',
#        x = '',
#        y = 'Average Survival') +
#   theme_classic() +
#   theme(axis.title = element_text(size = 12, face = "bold"),  # Larger axis labels
#         axis.text = element_text(size = 12),  # Readable tick labels
#         legend.title = element_text(size = 12, face = "bold"),
#         legend.text = element_text(size = 10),
#         plot.title = element_text(size = 12, face = "bold", hjust = 0)) # Centered title

ggsave('./figures/kelt_lgd_avg_survival_sex.png', fig_avg_s,
       width = 11, height = 8.5,
       dpi = 600, units = 'in', device = "png")

# save important results
save(mod_df, mod_dat, kelt.cjs.models, best_mod, best_res, phi_lgd,  file = paste0('./data/output/mod_dat_sy', max_sy, '_', Sys.Date(), '.rda'))

### END SCRIPT
