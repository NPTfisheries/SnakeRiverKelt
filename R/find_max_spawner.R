
# turn all obs prior to last spawner obs to spawner obs
# this often works, but sometimes errantly turns kelt observations to spawner observations. E.g., when a kelt briefly dips into a 
# downriver tributary (e.g., Walla Walla, John Day) before continuing downstream migration. Also, when a kelt has an errant time-stamp
# to a MRR site. At some MRR sites, sometime batch uploads occur and the min_det can be far later than when the fish likely 
# actually arrived at the MRR site.
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
