# Load model results and summarize.
# Author: Ryan N. Kinzer
# Date: 2024-05-31


library(tidyverse)

# load results-----
yr <- 2024
load(paste0('./data/cjs_model_data_',yr,'.Rda'))

# exact modeled data
mod_dat <- ch_mod_dat %>%
  filter(spawner_above == 1) %>%
  filter(!is.na(TRT_POPID)) %>%
  select(spawner_above, kelt_above, kelt_LGR, kelt_rs, group1 = spawn_year, group2 = TRT_POPID)

load(paste0('./data/mod_runs_above_',yr,'.Rda'))

#summary(mod.t.year.pop) #20 minutes at 2000 and 500
#plot(mod.t.year.pop)
#mod.t.year <- update(mod.t.year,n.iter = 1000)

mpg_df <- ch_mod_dat %>%
  select(MPG, POP_NAME, TRT_POPID) %>%
  mutate(pop = as.character(as.integer(factor(TRT_POPID)))) %>%
  distinct()

# model selection

dic.table <- data.frame("DIC" = sapply(model_ls$models, function(x){x$DIC}),
                        #"Rhat" = sapply(model_ls$models, function(x){x$Rhat}),
                        #"n.eff" = sapply(model_ls$models, function(x){x$n.eff}), 
                        "pD" = sapply(model_ls$models, function(x){x$pD}))

dic.table<-with(dic.table,dic.table[order(DIC),])
dic.table

# tmp <- mod.t$samples %>%
#   ggmcmc::ggs()
best_mod <- model_ls$models$mod.t.year.pop

mod_res <- best_mod$summary %>%
  as_tibble(rownames = 'rowId') %>%
  mutate(param = str_extract(rowId, '([^\\[]+)'))

param_df <- mod_res %>%
  filter(param %in% c('PHI', 'P')) %>%
  rowwise() %>%
  mutate(s = list(paste(strsplit(gsub('^.*\\[|]', '', rowId), NULL)[[1]], collapse=''))) %>%
  separate(s, into = c('pop', 'year','time')) %>%
  mutate(year = as.numeric(year) + 2009,#factor(year, labels = model_ls$data$grp1),         
         time = factor(time, labels = c('POP', 'LWG', 'Below')) 
         ) %>%
  left_join(mpg_df, by = 'pop')


# survivals
clear_kelt <- param_df %>%
  filter(param == 'PHI') %>%
  #filter(time == 2) %>%
  filter(time %in% c('POP')) %>%
  filter(MPG == 'Clearwater River') %>%
  ggplot(aes(x = fct_rev(as.factor(year)), y = `50%`, colour = TRT_POPID)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = .5)) +
  scale_color_brewer(palette = 'Dark2') +
  ylim(c(0,1)) +
  #facet_wrap(~MPG) +
  coord_flip() +
  labs(#title = 'Kelting Rates of Clearwater River Steelhead',
       #subtitle = 'Estimated probability of adult steelhead spawners leaving the population area as kelt.',
       x = 'Spawn Year',
       y = 'Proportion',
       colour = '') +
  theme_classic() +
  theme(legend.position = 'top')

clear_kelt

ggsave(paste0('./figures/clearwater_kelting_rates_',yr,'.png'), clear_kelt, width = 7, height = 7)


tmp <- param_df %>%
  group_by(param, time, MPG, TRT_POPID) %>%
  summarise(min_est = min(`50%`),
            mean_est = mean(`50%`),
            max_est = max(`50%`),)


clear_surv <- param_df %>%
  filter(param == 'PHI') %>%
  #filter(time == 2) %>%
  filter(time %in% c('LWG')) %>%
  filter(MPG == 'Clearwater River') %>%
  ggplot(aes(x = fct_rev(as.factor(year)), y = `50%`, colour = TRT_POPID)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = .5)) +
  scale_color_brewer(palette = 'Dark2') +
  ylim(c(0,1)) +
  facet_wrap(~MPG) +
  coord_flip() +
  labs(#title = 'Surival of Clearwater River Steelhead Kelt',
       #subtitle = 'Estimated probability of adult steelhead kelt leaving their spawning area as kelt \n and suriving to Lower Granite Dam.',
       x = 'Spawn Year',
       y = 'Survival',
       colour = '') +
  theme_classic() +
  theme(legend.position = 'top')

clear_surv

ggsave(paste0('./figures/clearwater_kelt_survival_',yr,'.png'), clear_surv, width = 7, height = 7)

# detection probs
clear_kelt_p <- param_df %>%
  filter(param == 'P') %>%
  filter(time %in% c('POP')) %>%
  filter(MPG == 'Clearwater River') %>%
  ggplot(aes(x = fct_rev(as.factor(year)), y = `50%`, colour = TRT_POPID)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = .5)) +
  scale_color_brewer(palette = 'Dark2') +
  ylim(c(0,1)) +
  #facet_wrap(~MPG) +
  coord_flip() +
  labs(#title = 'Detection of Clearwater River Steelhead Kelt',
       #subtitle = 'Estimated detection probability adult steelhead spawners \n leaving the population area as kelt.',
       x = 'Spawn Year',
       y = 'Proportion',
       colour = '') +
  theme_classic() +
  theme(legend.position = 'top')

clear_kelt_p

ggsave(paste0('./figures/clearwater_kelt_rate_detection_',yr,'.png'), clear_kelt_p, width = 7, height = 7)

clear_surv_p <- param_df %>%
  filter(param == 'P') %>%
  #filter(time == 2) %>%
  filter(time %in% c('LWG')) %>%
  filter(MPG == 'Clearwater River') %>%
  ggplot(aes(x = fct_rev(as.factor(year)), y = `50%`, colour = TRT_POPID)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = .5)) +
  scale_color_brewer(palette = 'Dark2') +
  ylim(c(0,1)) +
  #facet_wrap(~MPG) +
  coord_flip() +
  labs(#title = 'Detection of Clearwater River Steelhead Kelt',
       #subtitle = 'Estimated detection probability of steelhead kelt at Lower Granite Dam.',
       x = 'Spawn Year',
       y = 'Survival',
       colour = '') +
  theme_classic() +
  theme(legend.position = 'top')

clear_surv_p

ggsave(paste0('./figures/clearwater_kelt_survival_detection_',yr,'.png'), clear_surv_p, width = 7, height = 7)


# all survival to LWG
param_df %>%
  filter(param == 'PHI') %>%
  filter(time %in% c('LWG')) %>%
  ggplot(aes(x = TRT_POPID, y = `50%`, colour = `50%`)) +
  geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`), position = position_dodge(width = .5)) +
  #scale_color_brewer(palette = 'Dark2') +
  ylim(c(0,1)) +
  #facet_wrap(~year) +
  coord_flip() +
  labs(#title = 'Surival of Clearwater River Steelhead Kelt',
    #subtitle = 'Estimated probability of adult steelhead kelt leaving their spawning area as kelt \n and suriving to Lower Granite Dam.',
    x = 'Spawn Year',
    y = 'Survival',
    colour = '') +
  theme_classic() +
  theme(legend.position = 'top')




names(best_mod)


mod.sims.phi <- best_mod$sims.list$PHI %>%
  as_tibble() %>%
  gather(key = 'time', value = phi) %>%
  mutate(time = as.integer(gsub('V','',time))) %>%
  group_by(time) %>%
  mutate(iter = 1:n())

ggplot(tmp, aes(x = value, colour = Chain, fill = Chain)) +
  geom_freqpoly() +
  facet_wrap(~Parameter, scales = 'free')
