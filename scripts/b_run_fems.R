# Purpose: Examine proportions and abundance of B-run size and female individuals
#   in select Snake River populations. I.e., Where would we most likely encounter larger females?
#
# Author: Ryan N. Kinzer
# Date Created: July 23, 2025
#   Last Modified: 

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(readxl)
library(ggbreak)

b_pops = c("CRLMA-s/CRSFC-s",
           "CRLOC-s",
           "CRLOL-s",
           "CRSEL-s",
           #"CRSFC-s",
           "SFMAI-s",
           #"SFMAI-s/SFSEC-s",
           "SFSEC-s")

#--------------
# load data
pop_esc_df = read_xlsx(path = "../SnakeRiverFishStatus/output/syntheses/LGR_Steelhead_all_summaries_2025-07-11.xlsx",
                       sheet = "Pop_Tot_Esc")

pop_sex_df = read_xlsx(path = "../SnakeRiverFishStatus/output/syntheses/LGR_Steelhead_all_summaries_2025-07-11.xlsx",
                       sheet = "Pop_Sex_Props")

pop_size_df = read_xlsx(path = "../SnakeRiverFishStatus/output/syntheses/LGR_Steelhead_all_summaries_2025-07-11.xlsx",
                        sheet = "Pop_Size_Props")

# df to calculate b-run size female abundance
b_f_df = pop_esc_df %>%
  filter(popid %in% b_pops) %>%
  select(species,
         spawn_yr,
         mpg,
         popid,
         pop_sites,
         n_tags,
         median,
         median_exp) %>%
  left_join(pop_sex_df %>%
              filter(param == "p_fem") %>%
              select(species,
                     spawn_yr,
                     popid,
                     p_fem = median),
            by = join_by(species, spawn_yr, popid)) %>%
  left_join(pop_size_df %>%
              filter(param == "p_b") %>%
              select(species,
                     spawn_yr,
                     popid,
                     p_b = median),
            by = join_by(species, spawn_yr, popid)) %>%
  mutate(popid = recode(popid, "CRLMA-s/CRSFC-s" = "CRSFC-s")) %>%
  mutate(n_fem_b = median_exp * (p_fem * p_b))

# line graph by population
b_f_p = b_f_df %>%
  ggplot(aes(x = spawn_yr,
             y = n_fem_b,
             color = popid,
             group = popid)) +
  geom_line(na.rm = FALSE) +
  geom_point() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -60, hjust = 0, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  ) +
  scale_x_continuous(breaks = unique(b_f_df$spawn_yr)) +
  scale_y_break(c(750, 1250),
                scales = 0.25) +
  labs(x = "Spawn Year",
       y = "Number of B-run Size Females",
       color = "TRT Population")
b_f_p

# calculate mean n_fem_b per popid for spawn_yr >= 2018
pop_order = b_f_df %>%
  filter(spawn_yr >= 2018) %>%
  group_by(popid) %>%
  summarise(mean_n_fem_b = mean(n_fem_b, na.rm = TRUE)) %>%
  arrange(mean_n_fem_b) %>%
  pull(popid)

# dot plot by population
b_f_p2 = b_f_df %>%
  filter(spawn_yr >= 2018) %>%
  mutate(popid = factor(popid, levels = pop_order)) %>%
  ggplot(aes(x = popid, y = n_fem_b, color = as.factor(spawn_yr))) +
  geom_jitter(width = 0, size = 3, alpha = 0.7) +  
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "black") +  # mean dot
  theme_minimal() +
  labs(
    x = NULL,
    y = "Number of B-run Size Females",
    color = "Spawn Year"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  )
b_f_p2

# save dot plot
ggsave("./figures/b_run_fem_abund.pdf",
     plot = b_f_p2,
     width = 8,
     height = 6)
