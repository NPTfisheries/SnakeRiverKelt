# Purpose: Turn all observations prior to the last spawner observation to a spawner observation. The hope is to avoid errant downstream
#   observations (e.g., an adult moves to a hydrosystem site after release at LGR, then later moves back upriver to spawn).
#
# Author: Ryan Kinzer
#
# Notes: This function often works, but sometimes errantly turns kelt observations to spawner observations. For example, a kelt may
# briefly dip into a downriver tributary (e.g., Walla Walla, John Day) before continuing downstream migration as a kelt. Also, a kelt may
# have an errant time-stamp to a MRR site. At some MRR sites, sometimes batch uploads occur and the min_det date-time can be far later
# than when that fish actually arrived at the MRR iste
#
#' Find max spawner observation.
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
find_max_spawner = function(df) {
  tmp = df %>%
    filter(life_stage == "spawner") %>%
    group_by(tag_code) %>%
    slice(which.max(min_det)) %>%
    select(tag_code, max_spawner_det = min_det) %>%
    right_join(df) %>%
    mutate(tmp = case_when(
      min_det <= max_spawner_det ~ "spawner",
      TRUE ~ life_stage
    ))
  
  return(tmp)
}
