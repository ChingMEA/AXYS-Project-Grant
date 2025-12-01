# inter_lab_variability_helpers.r

# ----- Outline -----
# icc2_from_psych
# get_vc_row
# vc_props

# ---- icc2_from_psych ----
icc2_from_psych <- function(icc_obj) {
  # psych::ICC() returns a list; $results is a data.frame with rows of ICC types
  row <- icc_obj$results[icc_obj$results$type == "ICC2", ]
  
  list(
    ICC2      = as.numeric(row$ICC),
    ICC2_low  = as.numeric(row$`lower bound`),
    ICC2_high = as.numeric(row$`upper bound`),
    ICC2_p    = as.numeric(row$p)
  )
}

# ---- get_vc_row ----
# Pull % variance for a given grouping factor (e.g., "Donor", "Lab", "Residual") from a variance components data.frame produced by as.data.frame(VarCorr(mod))
get_vc_row <- function(vc_df, grp_name) {
  vc_df$prop_pct[vc_df$grp == grp_name]
}

# ---- vc_props ----
# Compute variance components (absolute and %) for a given lmer model
# Returns a tibble with Metric, Type (e.g., "Raw", "Standardized"), grouping factor (grp), variance (vcov), and percent of total variance
vc_props <- function(mod, metric, type) {
  vc <- as.data.frame(VarCorr(mod)) |>
    dplyr::select(grp, vcov)
  
  total <- sum(vc$vcov, na.rm = TRUE)
  
  vc |>
    dplyr::mutate(
      prop      = vcov / total,
      prop_pct  = 100 * prop,
      Metric    = metric,
      Type      = type
    ) |>
    dplyr::select(Metric, Type, grp, vcov, prop_pct)
}