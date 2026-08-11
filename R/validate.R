# Fail-fast configuration checks.
#
# These run once at load, before any modelling, because the joins they
# protect fail SILENTLY. A band named in the config but absent from the
# data produces NA waning at the left_join in apply_scenario(), which
# then propagates NA admissions all the way into the submission file
# without ever raising an error.


# Validate every age band referenced by the config against the bands
# actually present in the data.
#
# The check is deliberately asymmetric:
#
#   config -> data   A band named in the config but missing from the
#                    data is almost always a typo, and is a hard ERROR.
#
#   data -> config   A band in the data with no waning entry legitimately
#                    means "this programme does not reach this band". It
#                    defaults to 0 and is only REPORTED.
#
# Also enforces that the infant and adult programmes cover disjoint
# bands: the two are independent programmes, and a band claimed by both
# would have its coverage counted twice.
validate_age_groups <- function(cfg, epi) {

  data_bands <- sort(unique(epi$burden$age_group))

  # ---- config -> data: hard error ----
  referenced <- list(
    "infant_vaccination.waning_by_band"     = cfg$infant$waning_df$age_group,
    "infant_vaccination.age_bounds"         = cfg$infant$age_bounds$age_group,
    "adult_vaccination.eligible_age_groups" = cfg$adult$eligible_age_groups,
    "age_group_order"                       = unlist(cfg$age_group_order)
  )

  offenders <- Filter(length, lapply(referenced, setdiff, y = data_bands))

  if (length(offenders) > 0) {
    detail <- paste0("  ", names(offenders), ": ",
                     vapply(offenders,
                            function(x) paste0('"', x, '"', collapse = ", "),
                            character(1)),
                     collapse = "\n")
    stop("Age groups referenced in the config are not present in the data:\n",
         detail,
         "\nAge groups found in the data:\n  ",
         paste(data_bands, collapse = ", "),
         call. = FALSE)
  }

  # ---- programme overlap: hard error ----
  overlap <- intersect(cfg$infant$age_bounds$age_group,
                       cfg$adult$eligible_age_groups)
  if (length(overlap) > 0) {
    stop("The infant and adult programmes are independent and must cover ",
         "disjoint age groups.\n",
         "  Claimed by both: ", paste0('"', overlap, '"', collapse = ", "),
         "\n  infant_vaccination.age_bounds and ",
         "adult_vaccination.eligible_age_groups must not intersect.",
         call. = FALSE)
  }

  # ---- data -> config: report only ----
  no_waning <- setdiff(data_bands, cfg$infant$waning_df$age_group)
  if (length(no_waning) > 0) {
    message("Age groups with no infant waning entry (defaulting to 0): ",
            paste(no_waning, collapse = ", "))
  }

  invisible(TRUE)
}


# Validate adult-programme settings that are not covered by
# load_waning_curves() (which checks the ensemble file itself).
validate_adult_config <- function(cfg) {

  allowed <- c("zero", "hold_last")
  if (!isTRUE(cfg$adult$ve_beyond_curve %in% allowed)) {
    stop("adult_vaccination.ve_beyond_curve must be one of: ",
         paste(allowed, collapse = ", "),
         "\n  got: ", deparse(cfg$adult$ve_beyond_curve),
         call. = FALSE)
  }

  if (length(cfg$adult$eligible_age_groups) == 0) {
    message("adult_vaccination.eligible_age_groups is empty - ",
            "no adult programme will be applied.")
  }

  invisible(TRUE)
}


# Convenience wrapper: run every check.
validate_config <- function(cfg, epi) {
  validate_age_groups(cfg, epi)
  validate_adult_config(cfg)
  invisible(TRUE)
}
