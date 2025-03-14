# Create figures
# Author: Ryan N. Kinzer
# Date: 2024-06-04


library(tidyverse)
library(readxl)


# Population abundance----
# SFCLW Steelhead abundance compared to other populations.
pop_df <- read_xlsx(path = "../SnakeRiverFishStatus/output/syntheses/LGR_Steelhead_all_summaries_2024-04-19.xlsx",
                        sheet = "Pop_Tot_Esc") %>%
  filter(valid_est == 1) %>%
  group_by(TRT_POPID) %>%
  mutate(mu = mean(median,na.rm = TRUE),
         stdd = sd(median,na.rm = TRUE),
         z = (median - mu)/stdd,
         gm = exp(sum(log(median[median > 0]), na.rm = T) / length(median)))


pop_df %>%
  mutate(sfclw = ifelse(TRT_POPID == 'CRSFC-s', TRUE, FALSE)) %>%
  ggplot(aes(x = spawn_yr, 
             y = z, 
             group = TRT_POPID,
             colour = sfclw)) +
  #geom_ribbon(aes(ymin = lower95ci, ymax = upper95ci), alpha = .25) +
  geom_line() +
  #geom_hline(aes(yintercept = gm), linetype = 2) +
  #geom_hline(aes(yintercept = mat), color = 'red', linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = scales::breaks_pretty(3)) +
  scale_color_manual(values = c('grey70', 'black')) +
  # facet_wrap(~ TRT_POPID, 
  #            scales = "free_y", 
  #            ncol = 4, 
  #            labeller = label_wrap_gen(width = 19)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Spawn Year",
       y = "Population Escapement") +
  theme_bw() 

