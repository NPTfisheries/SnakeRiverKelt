# Create capture histories for steelhead including spawner, kelt, and repeat spawner observations.

# Purpose: A function to create capture histories for steelhead including spawner, kelt, and repeat spawner observations
#   using Ackerman's method.
#
# Author: Mike Ackerman
# Date Created: 03/17/2025
#   Date Last Modified: 03/17/2025
#   Modified By: Mike Ackerman
steelheadCapHist_MA <- function(x) {
  
  ch = x %>%
    summarise(
      release_lgr = if_else(any(max_det == lgr_max_det & life_stage == "spawner" & node == "LGR"), 1, 0),
      spawner_above = if_else(any(str_starts(rkm, "522") & rkm_total > 695 & min_det > lgr_max_det & life_stage == "spawner"), 1, 0),
      spawner_below = if_else(any(life_stage == "spawner" & grepl("GRS", path) & min_det > lgr_max_det & !site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON")), 1, 0),
      kelt_above = if_else(any(life_stage == "kelt"    & !grepl("GRS", path) & node != "LGR"), 1, 0),
      kelt_grs = if_else(any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      kelt_goa = if_else(any(site_code == "GOA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      kelt_lma = if_else(any(site_code == "LMA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      kelt_ihr = if_else(any(site_code == "IHR" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      kelt_mcn = if_else(any(site_code == "MCN" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      kelt_jda = if_else(any(site_code == "JDA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      kelt_tda = if_else(any(site_code == "TDA" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      kelt_bon = if_else(any(site_code == "BON" & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      #kelt_dwn = if_else(any(site_code %in% c("GOA", "LMA", "IHR", "MCN", "JDA", "TDA", "BON") & life_stage == "kelt" & min_det > lgr_max_det), 1, 0),
      rs_bon   = if_else(any(life_stage == "repeat spawner" & node == "BON"), 1, 0),
      rs_lgr   = if_else(any(life_stage == "repeat spawner" & node == "LGR"), 1, 0),
      rs_above = if_else(any(life_stage == "repeat spawner" & !grepl("GRS", path) & node != "LGR"), 1, 0),
      # when was the fish last observed moving upstream at lgr?
      lgr_max_det = unique(lgr_max_det),
      # if grs = 1, when was the fish last observed at grs as a kelt?
      grs_kelt_det = if (any(site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det)) {
        max(max_det[site_code == "GRS" & life_stage == "kelt" & min_det > lgr_max_det], na.rm = TRUE)
      } else { NA },
      .groups = "drop"
    )
  
  return(ch)
  
} # end steelheadCapHist_MA()
  