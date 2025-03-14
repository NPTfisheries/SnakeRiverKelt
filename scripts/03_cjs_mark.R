# Purpose: Run preliminary models on single populations/years using Mark to
# test other Bayesian model outputs for validation.
# Author: Ryan N. Kinzer
# Date: 2021-10-08

# load pkgs----
library(tidyverse)
library(RMark)

# load data----
yr <- 2024

load(paste0('./data/cjs_model_data_',yr,'.Rda'))

ch_data <- ch_mod_dat %>%
  filter(spawner_above == 1) %>%
  filter(!is.na(TRT_POPID)) %>%
  select(ch = ch_above, everything()) %>%
  mutate(pop = as.character(as.numeric(factor(TRT_POPID))))

#use ch_filtered at first for debugging
kelt.processed=process.data(ch_data, model="CJS", begin.time=1)#, group = 'spawn_year')
kelt.ddl=make.design.data(kelt.processed)
#null_model<-mark(kelt.processed,kelt.ddl)

Phi = list(formula=~1)
Phi.t = list(formula=~time)
Phi.year = list(formula=~spawn_year+time)
Phi.pop = list(formula=~pop+time)

p = list(formula=~1)
p.t = list(formula=~time)
p.year = list(formula=~spawn_year+time)
p.pop = list(formula=~pop+time)

kelt.model.list <- create.model.list("CJS")
kelt.results <- mark.wrapper(kelt.model.list,
                              data = kelt.processed, ddl = kelt.ddl, delete = TRUE)

model_summary <- kelt.results$model.table

best_mod <- kelt.results$Phi.pop.p.year

#mod.t <- mark(kelt.processed, kelt.ddl, model.parameters=list(Phi=Phi.t, p=p.t))
#mod.year.t <- mark(kelt.processed, kelt.ddl, model.parameters=list(Phi=Phi.year.t, p=p.year.t))
#mod.pop <- mark(kelt.processed, kelt.ddl, model.parameters=list(Phi=Phi.pop, p=p.pop))
#mod1 <- mark(kelt.processed, kelt.ddl, model.parameters=list(Phi=Phi.dot, p=p.dot)) # equal S and p across pops and time
#mod2 <- mark(kelt.processed, kelt.ddl, model.parameters=list(Phi=Phi.pop, p=p.pop)) # diff S and p across pops
#mod3 <- mark(kelt.processed, kelt.ddl, model.parameters=list(Phi=Phi.pop.t, p=p.pop.t), output = FALSE) # diff S and p across pops and time

#PIMS(mod3, "Phi", simplified = TRUE)
#results <- collect.models()
#results
#summary(mod1)
#summary(mod2)
#summary(mod3)
#https://sites.google.com/site/uwrandbayesianmodeling/schedule/7-capture-mark-recapture-modeling/7-3-open-dynamic-cmr-models/7-3-1-cjs-model-in-r

mod_ests <- best_mod$results$real %>%
  mutate(grp = rownames(.)) %>%
  mutate(param = str_split(grp, ' ',simplify = TRUE)[,1],
         TRT_POPID = str_split(grp, ' g', simplify = TRUE)[,2],
         TRT_POPID = str_split(TRT_POPID, ' ', simplify = TRUE)[,1],
         time = str_sub(grp, start = -1, end = -1)) %>%
  # mutate(TRT_POPID = c(rep(mod3$group.labels,each = 3), rep(mod3$group.labels,each = 3)),
  #        TRT_POPID = gsub('TRT_POPID', '', TRT_POPID)) %>%
  select(-grp) %>%
  select(TRT_POPID, param, time, estimate:ucl)


mod_ests %>%
  #filter(time == 1) %>%
  filter(param == 'Phi') %>%
  ggplot(aes(x = fct_rev(time), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
  geom_point() +
  coord_flip() +
  labs(x = '',
       y = 'Survival Probability')



sumod3_ests %>%
  filter(time == 2) %>%
  filter(param == 'Phi') %>%
  ggplot(aes(x = fct_reorder(st_POP_NAME, estimate), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
  geom_point() +
  coord_flip() +
  labs(x = '',
       y = 'Survival Probability') +
  theme_rk()


mod3_ests %>%
  filter(time == 2) %>%
  filter(param == 'p') %>%
  ggplot(aes(x = fct_reorder(st_POP_NAME, estimate), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
  geom_point() +
  coord_flip() +
  labs(x = '',
       y = 'Detection Probability') +
  theme_rk()

mod3_ests %>%
  filter(time == 3) %>%
  filter(param == 'p') %>%
  ggplot(aes(x = fct_reorder(st_POP_NAME, estimate), y = estimate)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
  geom_point() +
  coord_flip() +
  labs(x = '',
       y = 'Detection Probability') +
