# load packages
library(tidyverse)
library(marked)

# set arguments
max_sy <- 2024

# load data
load(paste0('./data/input/cjs_model_data_sy', max_sy,'.Rda'))

df = ch_mod_dat %>%
  # keep just records with spawning observations and a final spawning population
  filter(spawner_above == 1 & !is.na(popid)) %>%
  # recode some spawners to populations, if desired
  mutate(popid = recode(popid, "CRLMA-s/CRSFC-s" = "CRSFC-s")) %>%
  # remove some observations with an unknown final population (e.g., SFG, USE, USI)
  filter(!grepl("/", popid))

#df <- ch_mod_dat[ch_mod_dat$spawner_above == 1 & !is.na(ch_mod_dat$popid) & !grepl("/", ch_mod_dat$popid),]

df %>% 
  group_by(spawn_yr) %>% 
  summarise(spawner_above = sum(spawner_above),
            kelt_above = sum(kelt_above),
            kelt_grs = sum(kelt_grs),
            kelt_below = sum(kelt_below))


mod_dat <- data.frame(
  ch = paste0(df$spawner_above, df$kelt_grs, df$kelt_below),
  freq = rep(1, dim(df)[1]),
  mpg = df$mpg,
  pop = df$popid,
  year = df$spawn_yr
)

mod_dat$mpg <- as.factor(mod_dat$mpg)
mod_dat$pop <- as.factor(mod_dat$pop)
mod_dat$year <- as.factor(mod_dat$year)

# create function to run all models
fit.kelt.cjs.models <- function(){
  
  # Apparent survival (Phi) formula
  Phi.mpg.time <- list(formula=~mpg*time)
  Phi.pop.time <- list(formula=~pop*time)  
  Phi.year.time <- list(formula=~year*time)
  Phi.year.mpg <- list(formula=~year+mpg+time)
  Phi.year.pop <- list(formula=~year+pop+time)
  Phi.mpg <- list(formula=~mpg + time)
  Phi.pop <- list(formula=~pop + time)
  Phi.year <- list(formula=~year + time)
  Phi.time <- list(formula=~time) 
  #Phi.dot <- list(formula=~1) # constant survival
  
  # Detection probability (p) formula
  #p.pop.time <- list(formula=~pop*time) 
  p.year.time <- list(formula=~year*time)
  #p.pop <- list(formula=~pop + time)  
  p.year <- list(formula=~year + time)
  p.time <- list(formula=~time) # one discrete estimate of p per capture event
  #p.dot <- list(formula=~1) # constant detection
  
  # Construct all combinations and put into one model table
  cml <- create.model.list(c("Phi","p")) # makes all possibile combinations of those parameter formulas
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

kelt.proc <- process.data(mod_dat, begin.time = 1, model = 'CJS', groups = c("mpg", "pop", "year"))
kelt.ddl <- make.design.data(kelt.proc)

kelt.cjs.models <- fit.kelt.cjs.models()

kelt.cjs.models
best_mod <- kelt.cjs.models[[24]] #$Phi.pop.p.pop.time

best_res <- best_mod$results$reals %>%
  bind_rows(.id = 'param')

# marginize over the year covariate to calcualte an average survival for each pop


source('./R/theme_rk.R')

fig_s <- best_res %>%
  filter(time == 1) %>%
  filter(param == 'Phi') %>%
  #filter(estimate < .99) %>%
  ggplot(aes(x = fct_reorder(pop,estimate), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
  geom_point() +
  coord_flip() +
  facet_wrap(~year, ncol= 5) +
  labs(title = 'Snake River Basin Steelhead Kelt',
       subtitle = 'Survival of adult steelhead from Snake River basin spawning areas to the kelt life-stage at Lower Granite Dam.',
       x = '',
       y = 'P(Survival | Observed Spawner)') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),  # Larger axis labels
        axis.text = element_text(size = 8),  # Readable tick labels
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0)) # Centered title

fig_s

ggsave('./figures/kelt_lgd_survival.png', fig_s,
       width = 11, height = 8.5,
       dpi = 600, units = 'in', device = "png")

fig_d <- best_res %>%
  filter(time == 2) %>%
  filter(param == 'p') %>%
  ggplot(aes(x = fct_rev(year), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
  geom_point() +
  coord_flip() +
  labs(title = 'Snake River Basin Steelhead Kelt',
       subtitle = 'Detection probability of steelhead kelt at Lower Granite Dam after being observed in Snake River basin spawning areas.',
       x = 'Spawn Year',
       y = 'P(Detection | Kelt at LGD)') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),  # Larger axis labels
        axis.text = element_text(size = 8),  # Readable tick labels
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0)) # Centered title

fig_d

ggsave('./figures/kelt_lgd_detection.png', fig_d,
       width = 11, height = 8.5,
       dpi = 600, units = 'in', device = "png")

# pop average - uses an inverse weightd mean to account for precision around each individual estimate
phi_lgd <- best_res %>%
  filter(time == 1) %>%
  filter(param == 'Phi') %>%
  group_by(pop, time) %>%
  summarize(
    avg_Phi = sum(estimate / (se^2), na.rm = TRUE) / sum(1 / (se^2), na.rm = TRUE),
    se_Phi = sqrt(1 / sum(1 / (se^2), na.rm = TRUE)),  # Weighted SE
    lower_CI = pmax(avg_Phi - 1.96 * se_Phi, 0),  # 95% CI lower bound
    upper_CI = pmin(avg_Phi + 1.96 * se_Phi, 1),  # 95% CI upper bound
    .groups = "drop"
  )

fig_avg_s <- phi_lgd %>%
  ggplot(aes(x = fct_reorder(pop,avg_Phi), y = avg_Phi)) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by = .25)) +
  geom_point() +
  coord_flip() +
  labs(title = 'Snake River Basin Steelhead Kelt',
       subtitle = 'Average survival of adult steelhead from Snake River population spawning areas to the kelt life-stage at Lower Granite Dam.',
       x = '',
       y = 'Average Survival') +
  theme_classic() +
  theme(axis.title = element_text(size = 12, face = "bold"),  # Larger axis labels
        axis.text = element_text(size = 12),  # Readable tick labels
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0)) # Centered title

fig_avg_s

ggsave('./figures/kelt_lgd_avg_survival.png', fig_avg_s,
       width = 11, height = 8.5,
       dpi = 600, units = 'in', device = "png")

save(df, mod_dat, kelt.cjs.models, best_mod, best_res, phi_lgd,  file =paste0('./data/output/mod_dat_sy', max_sy,'.Rda'))
