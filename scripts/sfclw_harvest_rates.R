# Purpose: Prepare a figure of catch/release numbers
#
# Author: Ryan N. Kinzer
# Date Created: 2025-04-23
#   Date Last Modified:
#   Modified By:

# load packages
library(tidyverse)

# set arguments
sy = 2025

sheets <- readxl::excel_sheets(paste0('./data/input/Steelhead_Harvest_Summary_SY',sy,'.xlsx'))

week_df <- tibble(week = sheets,
                  week_num = 1:length(sheets))

dat <- map_dfr(sheets,
               .f = ~readxl::read_xlsx(paste0('./data/input/Steelhead_Harvest_Summary_SY',sy,'.xlsx'),
                                       sheet = .,
                                       col_types = 'text',
                                       skip = 1),
               .id = 'week_num'
               )

names(dat) <- gsub(' ','_', tolower(names(dat)))


dat <- dat %>%
  filter(!is.na(section)) %>%
  mutate(week_num = as.integer(week_num),
         across(anglers_interviewed:total_catch, ~as.integer(.x))) %>%
  left_join(week_df, by = 'week_num') %>%
  mutate(p_kept = steelhead_kept/total_catch)


dat %>%
  ggplot(aes(x = week_num, y = p_kept)) +
  geom_line() +
  facet_wrap(~location)

sf_clw <- dat %>%
  filter(section == 7)

sf_clw %>%
  summarise(caught = sum(total_catch, na.rm = TRUE),
            released = sum(steelhead_released, na.rm = TRUE),
            kept = sum(steelhead_kept, na.rm = TRUE))

sf_clw %>%
ggplot(aes(x = week_num)) +
  geom_col(aes(y = total_catch), fill = "steelblue", alpha = 0.6) +  # total fish caught
  geom_line(aes(y = p_kept * max(total_catch, na.rm = TRUE), group = 1), color = "darkred", size = 1.2) +  # scale proportion to match y-axis
  geom_point(aes(y = p_kept * max(total_catch, na.rm = TRUE)), color = "darkred", size = 3) +
  scale_x_continuous(labels = sf_clw$week,
                     breaks = sf_clw$week_num) +
  scale_y_continuous(
    name = "Total Fish Caught",
    expand = c(0,0),
    sec.axis = sec_axis(
      transform = ~ (. / max(sf_clw$total_catch, na.rm = TRUE))*100,
      name = "Percent Kept")
  ) +
  labs(x = "Week", title = "SF Clearwater River",
  subtitle = "Weekly Fish Catch and Proportion Kept") +
  theme_minimal() +
  theme(plot.background = element_rect(fill = 'white'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave('./figures/SFCLW_harvest_rates.png', width = 9, height = 7, dpi = 300)


sf_clw %>%
  select(week_num, week, contains('steelhead')) %>%
  pivot_longer(-c(week_num, week), names_to = 'key', values_to = 'value') %>%
  ggplot(aes(x = week_num)) +
  geom_col(aes(y = value, fill = key)) +
  scale_x_continuous(labels = sf_clw$week,
                     breaks = sf_clw$week_num) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c('darkred','steelblue'), labels = c('Kept', 'Released')) +
  labs(x = "Week",
       y = 'Total Fish Caught',
       fill = '',
       title = "SF Clearwater River",
       subtitle = "Weekly Fish Kept and Released") +
  theme_minimal() +
  theme(legend.position = 'top',
        plot.background = element_rect(fill = 'white'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave('./figures/SFCLW_catch.png', width = 9, height = 7, dpi = 300)
