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
validate_age_groups <- function(cfg, epi, quiet = FALSE) {

  data_bands <- sort(unique(epi$burden$age_group))

  # ---- config -> data: hard error ----
  referenced <- list(
    "infant_vaccination.waning_by_band"     = cfg$infant$waning_df$age_group,
    "infant_vaccination.age_bounds"         = cfg$infant$age_bounds$age_group,
    "adult_vaccination.eligible_age_groups" = cfg$adult$eligible_age_groups,
    "scenarios[].adult_age_groups"          = unlist(cfg$scenarios_df$adult_age_groups),
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

  # ---- data -> config: report only ----
  no_waning <- setdiff(data_bands, cfg$infant$waning_df$age_group)
  if (length(no_waning) > 0 && !quiet) {
    message("Age groups with no infant waning entry (defaulting to 0): ",
            paste(no_waning, collapse = ", "))
  }

  invisible(TRUE)
}


# Check every configured country resolves in all input files, before any
# modelling starts - so a typo fails immediately rather than 20 countries
# into a long run.
validate_countries <- function(cfg, raw) {

  cs <- cfg$countries_df

  if (nrow(cs) == 0) {
    stop("No countries configured. Add at least one entry under `countries`.",
         call. = FALSE)
  }
  if (anyDuplicated(cs$iso2) > 0) {
    dup <- unique(cs$iso2[duplicated(cs$iso2)])
    stop("Duplicate country iso2 code(s) in config: ",
         paste(dup, collapse = ", "), call. = FALSE)
  }
  if (any(!nzchar(cs$name)) || any(!nzchar(cs$iso2))) {
    stop("Every entry under `countries` needs a non-empty `name` and `iso2`.",
         call. = FALSE)
  }

  # target-data is keyed by full country name, auxiliary-data by ISO2
  for (k in c("admissions", "burden")) {
    avail <- unique(raw[[k]]$country)
    for (nm in cs$name) check_country_present(nm, avail, raw$labels[[k]])
  }
  for (k in c("births", "population")) {
    avail <- unique(raw[[k]]$country)
    for (cc in cs$iso2) check_country_present(cc, avail, raw$labels[[k]])
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

  # ---- campaigns ----
  camps <- cfg$adult$campaigns

  for (i in seq_along(camps)) {
    cp <- camps[[i]]
    if (is.na(cp$start) || is.na(cp$end)) {
      stop("adult_vaccination.campaigns[[", i, "]] has an unparseable ",
           "start or end date.", call. = FALSE)
    }
    if (cp$end < cp$start) {
      stop("adult_vaccination.campaigns[[", i, "]] ends before it starts: ",
           format(cp$start), " to ", format(cp$end), call. = FALSE)
    }
  }

  if (length(camps) > 0) {
    shares <- vapply(camps, function(c) c$share, numeric(1))
    if (any(is.na(shares))) {
      stop("adult_vaccination.campaigns: `share` must be given on every ",
           "campaign or omitted from all of them (for an equal split).",
           call. = FALSE)
    }
    if (abs(sum(shares) - 1) > 1e-8) {
      stop("adult_vaccination.campaigns: `share` values must sum to 1, ",
           "got ", format(sum(shares)),
           "\n  shares: ", paste(shares, collapse = ", "),
           call. = FALSE)
    }

    # One-off cumulative coverage assumes campaigns recruit distinct
    # people; overlapping windows make that ambiguous.
    if (length(camps) > 1) {
      ord <- order(vapply(camps, function(c) as.numeric(c$start), numeric(1)))
      s   <- camps[ord]
      for (i in seq_len(length(s) - 1)) {
        if (s[[i + 1]]$start <= s[[i]]$end) {
          warning("adult_vaccination.campaigns: windows ", i, " and ", i + 1,
                  " overlap. Coverage is one-off and cumulative, so ",
                  "overlapping campaigns are ambiguous.", call. = FALSE)
        }
      }
    }
  }

  if (!isTRUE(cfg$adult$baseline_coverage >= 0 &&
              cfg$adult$baseline_coverage <= 1)) {
    stop("adult_vaccination.baseline_coverage must lie in [0, 1]; got ",
         deparse(cfg$adult$baseline_coverage), call. = FALSE)
  }

  invisible(TRUE)
}


# Validate the scenario table.
validate_scenarios <- function(cfg) {

  sc <- cfg$scenarios_df

  if (nrow(sc) == 0) {
    stop("No scenarios configured. Add at least one entry under `scenarios`.",
         call. = FALSE)
  }
  if (anyDuplicated(sc$id) > 0) {
    stop("Duplicate scenario id(s): ",
         paste(unique(sc$id[duplicated(sc$id)]), collapse = ", "),
         call. = FALSE)
  }
  if (any(!nzchar(sc$id))) {
    stop("Every scenario needs a non-empty `id`.", call. = FALSE)
  }

  for (col in c("infant_uptake", "adult_coverage")) {
    v   <- sc[[col]]
    bad <- which(is.na(v) | v < 0 | v > 1)
    if (length(bad) > 0) {
      stop("scenarios: `", col, "` must lie in [0, 1].\n  Offending: ",
           paste0(sc$id[bad], " = ", v[bad], collapse = ", "),
           call. = FALSE)
    }
  }

  # A scenario asking for adult coverage with no bands to apply it to
  # would silently have no effect.
  empty_bands <- vapply(sc$adult_age_groups, length, integer(1)) == 0
  offenders   <- sc$id[sc$adult_coverage > 0 & empty_bands]
  if (length(offenders) > 0) {
    stop("Scenario(s) set adult_coverage > 0 but target no age groups, so ",
         "the coverage would have no effect: ",
         paste0('"', offenders, '"', collapse = ", "),
         "\n  Give them an `adult_age_groups` list, or populate ",
         "adult_vaccination.eligible_age_groups.", call. = FALSE)
  }

  invisible(TRUE)
}


# Run the checks that do not depend on any one country's data. Called
# once, before the country loop.
# Every eligible adult band must appear in the population file, for every
# configured country. A missing band would silently shrink the dose
# denominator rather than erroring - the sum would just be over fewer
# bands, and the dose count would come out too low with no warning.
validate_population_coverage <- function(cfg, raw) {

  # The union across the default and every scenario override - any of
  # them could end up as a dose denominator.
  eligible <- all_adult_bands(cfg)
  if (length(eligible) == 0) return(invisible(TRUE))

  for (cc in cfg$countries_df$iso2) {
    have <- raw$population %>%
      filter(country == cc) %>%
      pull(age_group) %>%
      unique()
    missing <- setdiff(eligible, have)
    if (length(missing) > 0) {
      stop("Population data for \"", cc, "\" is missing age band(s) needed ",
           "as the adult dose denominator: ",
           paste0('"', missing, '"', collapse = ", "),
           "\n  Present for that country: ", paste(sort(have), collapse = ", "),
           "\n  file: ", raw$labels$population,
           call. = FALSE)
    }
  }
  invisible(TRUE)
}


# The two programmes must cover disjoint age bands.
#
# This is what makes bind_rows() of the two protection tables safe: a
# band claimed by both would get two rows per key, and the joins in
# apply_scenario() would multiply admissions rather than fail. It is
# checked against the UNION of every adult band referenced anywhere, so a
# per-scenario override cannot smuggle in a band the infant programme
# already claims.
#
# Config-only - no data needed - so it runs in the global pre-flight
# rather than per country.
validate_programme_disjoint <- function(cfg) {

  overlap <- intersect(cfg$infant$age_bounds$age_group, all_adult_bands(cfg))
  if (length(overlap) > 0) {
    claimed_by <- vapply(overlap, function(b) {
      where <- character(0)
      if (b %in% cfg$adult$eligible_age_groups) {
        where <- c(where, "adult_vaccination.eligible_age_groups")
      }
      hits <- cfg$scenarios_df$id[vapply(cfg$scenarios_df$adult_age_groups,
                                         function(g) b %in% g, logical(1))]
      if (length(hits) > 0) {
        where <- c(where, paste0("scenario(s) ", paste(hits, collapse = ", ")))
      }
      paste0('"', b, '" (', paste(where, collapse = "; "), ")")
    }, character(1))

    stop("The infant and adult programmes are independent and must cover ",
         "disjoint age groups.\n  Claimed by both:\n    ",
         paste(claimed_by, collapse = "\n    "),
         "\n  These bands are already in infant_vaccination.age_bounds.",
         call. = FALSE)
  }
  invisible(TRUE)
}


# data_start must not exclude the data the sampler needs.
#
# The burden file gives ONE season-long window per band, so the filter
# `burden_start_date >= data_start` is all-or-nothing: a data_start even
# one day after the window start empties the whole age-stratified table,
# leaving the sampler with no column margins and nothing to draw. That
# surfaced only as "no baseline samples produced", deep inside the
# country loop, so it is checked here instead.
validate_data_start <- function(cfg, raw) {

  ds <- ymd(cfg$data_start)
  if (is.na(ds)) {
    stop("data_start is not a valid date: ", deparse(cfg$data_start),
         call. = FALSE)
  }

  fmt_b <- cfg$input_date_formats$burden_agegroups
  fmt_w <- cfg$input_date_formats$weekly_counts

  for (nm in cfg$countries_df$name) {

    bur <- raw$burden     %>% filter(country == nm)
    adm <- raw$admissions %>% filter(country == nm)
    if (nrow(bur) == 0 || nrow(adm) == 0) next   # validate_countries reports this

    b_start <- min(parse_dates_strict(bur$start_date, fmt_b,
                                      "hospitalburden_agegroups.csv$start_date"))
    b_end   <- max(parse_dates_strict(bur$end_date, fmt_b,
                                      "hospitalburden_agegroups.csv$end_date"))
    a_rng   <- range(parse_dates_strict(adm$target_end_date, fmt_w,
                                        "hospitaladmissions.csv$target_end_date"))

    if (ds > b_start || ds > a_rng[2]) {
      stop("data_start (", format(ds), ") excludes the data needed for \"",
           nm, "\".",
           "\n  burden window:     ", format(b_start), " .. ", format(b_end),
           "\n  weekly admissions: ", format(a_rng[1]), " .. ", format(a_rng[2]),
           "\n",
           "\n  The burden file has a single season-long window per age band,",
           "\n  so any data_start after ", format(b_start), " drops the whole",
           "\n  age-stratified table and the sampler has nothing to draw.",
           "\n  Set data_start to ", format(b_start), " or earlier.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}


validate_global_config <- function(cfg, raw) {
  validate_countries(cfg, raw)
  validate_data_start(cfg, raw)
  validate_adult_config(cfg)
  validate_scenarios(cfg)
  validate_programme_disjoint(cfg)
  validate_population_coverage(cfg, raw)
  invisible(TRUE)
}
