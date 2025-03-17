# Create capture histories for steelhead including spawner, kelt, and repeat spawner observations.
#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
steelheadCapHist <- function(x) {
  
  # single pathway spawner observations and all kelt/repeat spawners
  ch <- x %>%
    #filter(life_stage != 'spawner') %>%
    mutate(obs_loc = case_when(
      life_stage == "spawner" & node == "LGR"                              ~ "release_lgr",   
      life_stage == "spawner" & !grepl("GRS", path) & node != "LGR"        ~ "spawner_above",
      life_stage == "spawner" &  grepl("GRS", path) & 
        !site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON") ~ "spawner_below",
      #life_stage == "spawner" &  grepl("GRS", path)                        ~ "spawner_below",
      life_stage == "kelt"    & !grepl("GRS", path) & node != "LGR"        ~ "kelt_above",
      life_stage == "kelt" & node %in% c("GRS")                            ~ "kelt_grs", # ma changed
      #life_stage == "kelt" & node %in% c("LGR", "GRS")                     ~ "kelt_grs",
      life_stage == "kelt" & node == "GOA"                                 ~ "kelt_goa",
      life_stage == "kelt" & node == "LMA"                                 ~ "kelt_lma",
      life_stage == "kelt" & node == "IHR"                                 ~ "kelt_ihr",
      life_stage == "kelt" & node == "MCN"                                 ~ "kelt_mcn",
      life_stage == "kelt" & node == "JDA"                                 ~ "kelt_jda",
      life_stage == "kelt" & node == "TDA"                                 ~ "kelt_tda",
      life_stage == "kelt" & node == "BON"                                 ~ "kelt_bon", 
      life_stage == "repeat spawner" & node == "BON"                       ~ "rs_bon",
      life_stage == "repeat spawner" & node == "LGR"                       ~ "rs_lgr",
      life_stage == "repeat spawner" & !grepl("GRS", path) & node != "LGR" ~ "rs_above",
      TRUE ~ "other")) %>%
    mutate(obs_loc = factor(obs_loc, levels = c('release_lgr', 'spawner_above', 'spawner_below', 'kelt_above', 'kelt_grs', 'kelt_goa', 'kelt_lma', 'kelt_ihr', 'kelt_mcn', 'kelt_jda', 'kelt_tda', 'kelt_bon', 'rs_bon', 'rs_lgr', 'rs_above', 'other'))) %>%
    arrange(obs_loc) %>%
    group_by(tag_code, obs_loc) %>%
    summarise(obs = 1) %>%
    pivot_wider(names_from = "obs_loc", values_from = "obs", values_fill = 0)
  
  col_order <- c("tag_code", "release_lgr", "spawner_above", "spawner_below", "kelt_above", "kelt_grs", "kelt_goa", "kelt_lma", "kelt_ihr", "kelt_mcn", "kelt_jda", "kelt_tda", "kelt_bon", "rs_bon", "rs_lgr", "rs_above", "other")
  
  if(!(all(col_order %in% colnames(ch)))){
    missed <- setdiff(col_order, colnames(ch))
    ch[missed] <- rep(0,nrow(ch))
  }
  
  ch <- ch %>%
    select(na.omit(match(col_order, names(.))))
  
  return(ch)
  
}