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
    "catchup_vaccination.age_bounds"        = cfg$catchup$age_bounds$age_group,
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

  for (col in c("adult_coverage", "catchup_coverage")) {
    v   <- sc[[col]]
    bad <- which(is.na(v) | v < 0 | v > 1)
    if (length(bad) > 0) {
      stop("scenarios: `", col, "` must lie in [0, 1].\n  Offending: ",
           paste0(sc$id[bad], " = ", v[bad], collapse = ", "),
           call. = FALSE)
    }
  }

  # infant_uptake is a VECTOR per scenario - one value per vaccination
  # window - so every element has to be checked, not just the first.
  bad <- which(vapply(sc$infant_uptake,
                      function(v) any(is.na(v) | v < 0 | v > 1), logical(1)))
  if (length(bad) > 0) {
    stop("scenarios: `infant_uptake` must lie in [0, 1] for every ",
         "vaccination window.\n  Offending: ",
         paste0(sc$id[bad], " = [",
                vapply(sc$infant_uptake[bad], paste, character(1),
                       collapse = ", "), "]", collapse = "; "),
         call. = FALSE)
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


# The infant and adult programmes must cover disjoint age bands.
#
# They are different vaccines given to different people, so a band
# claimed by both is a configuration error rather than a modelling
# choice. It is checked against the UNION of every adult band referenced
# anywhere, so a per-scenario override cannot smuggle in a band the
# infant programme already claims.
#
# The CATCH-UP programme is deliberately exempt. Its cohort ages through
# the infant bands, so it necessarily shares them with the birth-dose
# programme - that overlap is the point, not a mistake. What used to make
# overlap unsafe was bind_rows() giving two rows per key and the joins in
# apply_scenario() multiplying admissions; combine_protection() now sums
# the programmes within a band instead, and errors if their combined
# coverage exceeds 1.
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

    b_starts <- sort(unique(parse_dates_strict(
      bur$start_date, fmt_b, "hospitalburden_agegroups.csv$start_date")))
    b_end    <- max(parse_dates_strict(
      bur$end_date, fmt_b, "hospitalburden_agegroups.csv$end_date"))
    a_dates  <- sort(unique(parse_dates_strict(
      adm$target_end_date, fmt_w, "hospitaladmissions.csv$target_end_date")))

    # What matters is that BOTH tables still have rows after the filter -
    # not where the first window happens to start. A burden file with one
    # season-long window per band is emptied by any later data_start; one
    # with 4-weekly windows is not, and there a mid-series start is a
    # perfectly good way to select which seasons to model.
    kept_b <- b_starts[b_starts >= ds]
    kept_a <- a_dates[a_dates >= ds]

    if (length(kept_b) == 0 || length(kept_a) == 0) {

      why <- if (length(b_starts) == 1) {
        paste0("\n  The burden file has a SINGLE window per age band, so any",
               "\n  data_start after ", format(b_starts), " drops the whole",
               "\n  age-stratified table and the sampler has nothing to draw.",
               "\n  Set data_start to ", format(b_starts), " or earlier.")
      } else {
        paste0("\n  The burden file has ", length(b_starts), " windows per age ",
               "band. None begins on or", "\n  after data_start",
               if (length(kept_a) == 0)
                 ", and no weekly row survives it either" else "", ".",
               "\n  Set data_start to ", format(max(b_starts)), " or earlier.")
      }

      stop("data_start (", format(ds), ") excludes the data needed for \"",
           nm, "\".",
           "\n  burden windows:    ", format(min(b_starts)), " .. ",
           format(b_end), "  (", length(b_starts), " window(s), ",
           length(kept_b), " kept)",
           "\n  weekly admissions: ", format(min(a_dates)), " .. ",
           format(max(a_dates)), "  (", length(a_dates), " week(s), ",
           length(kept_a), " kept)",
           why, call. = FALSE)
    }
  }
  invisible(TRUE)
}


# Every adult campaign must overlap the modelled weekly grid.
#
# Campaign windows are matched against week-ending dates, so a window
# that lands on any other weekday matches nothing and silently delivers
# zero coverage - the run then completes normally and emits a full
# submission of zeros. Checked here so it fails before any modelling
# rather than surfacing later as an unrelated-looking error.
validate_campaign_windows <- function(cfg, raw) {

  weeks <- sort(unique(parse_dates_strict(
    raw$admissions$target_end_date, cfg$input_date_formats$weekly_counts,
    "hospitaladmissions.csv$target_end_date")))
  weeks <- weeks[weeks >= ymd(cfg$data_start)]
  if (length(weeks) == 0) return(invisible(TRUE))   # validate_data_start reports

  # Both campaign programmes are checked. A programme no scenario uses is
  # skipped: its window is then genuinely irrelevant.
  programmes <- list(
    list(key   = "adult_vaccination",
         camps = cfg$adult$campaigns,
         used  = any(cfg$scenarios_df$adult_coverage > 0)),
    list(key   = "catchup_vaccination",
         camps = cfg$catchup$campaigns,
         used  = any(cfg$scenarios_df$catchup_coverage > 0))
  )

  for (pr in programmes) {

    if (length(pr$camps) == 0 || !pr$used) next

    for (i in seq_along(pr$camps)) {
      cp <- pr$camps[[i]]
      if (!any(weeks >= cp$start & weeks <= cp$end)) {
        near <- weeks[order(abs(as.numeric(weeks - cp$start)))][1:2]
        stop(pr$key, ".campaigns[[", i, "]] (", format(cp$start),
             " to ", format(cp$end), ") matches no modelled week, so every ",
             "scenario would deliver zero coverage.",
             "\n  Campaign windows are matched against week-ending dates, ",
             "which are always ", weekdays(weeks[1]), "s.",
             "\n  ", format(cp$start), " is a ", weekdays(cp$start), ".",
             "\n  Nearest valid dates: ", paste(format(sort(near)), collapse = ", "),
             "\n  Widen the window, or move it onto a week-ending date.",
             call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}


# The horizon anchor must sit on the modelled weekly grid.
#
# horizon is (week - anchor) / 7 truncated toward zero by as.integer().
# When the anchor is OFF the grid, that truncation rounds negative
# differences up and positive ones down, so the week before the anchor
# and the week after it both land on horizon 0. Two weeks then share a
# horizon, and the submission carries a full week of duplicate rows for
# every scenario, draw and pop_group.
#
# Nothing about the data is wrong when that happens, so the failure
# surfaces far downstream as "duplicate rows" and an off-by-one grid
# count - which says nothing about the anchor. Hence this check.
validate_horizon_anchor <- function(cfg, raw) {

  weeks <- sort(unique(parse_dates_strict(
    raw$admissions$target_end_date, cfg$input_date_formats$weekly_counts,
    "hospitaladmissions.csv$target_end_date")))
  weeks <- weeks[weeks >= ymd(cfg$data_start)]
  if (length(weeks) == 0) return(invisible(TRUE))   # validate_data_start reports

  a <- cfg$anchor
  if (length(a) != 1 || is.na(a)) {
    stop("submission_horizon_anchor is missing or unparseable.", call. = FALSE)
  }

  if (any(as.numeric(weeks - a) %% 7 != 0)) {
    stop("submission_horizon_anchor (", format(a), ", a ", weekdays(a),
         ") is not on the modelled weekly grid, whose weeks end on ",
         weekdays(weeks[1]), "s.",
         "\n  horizon is (week - anchor) / 7 truncated toward zero, so an ",
         "off-grid anchor makes the week before it and the week after it ",
         "collide on horizon 0.",
         "\n  That surfaces later as duplicate submission rows rather than ",
         "as an anchor problem.",
         "\n  The first modelled week is ", format(min(weeks)), ".",
         call. = FALSE)
  }

  if (a != min(weeks)) {
    n <- as.integer((min(weeks) - a) / 7)
    stop("submission_horizon_anchor (", format(a), ") is not the first ",
         "modelled week (", format(min(weeks)), "), so horizon would ",
         if (n > 0) paste0("start at ", n) else paste0("start at ", n,
           " and run negative"),
         " rather than at 0.",
         "\n  Set it to ", format(min(weeks)), ", or move data_start.",
         call. = FALSE)
  }

  invisible(TRUE)
}


# Infant vaccination windows must be well-formed and non-overlapping.
#
# Overlap matters now that uptake is per window: a birth day inside two
# windows has two candidate uptakes, and both the coverage sum and the
# dose attribution would quietly pick one. Non-overlapping windows make
# the question unanswerable-by-construction instead.
validate_infant_windows <- function(cfg) {

  st <- cfg$infant$vacc_start
  en <- cfg$infant$vacc_end

  if (length(st) != length(en)) {
    stop("infant_vaccination.windows: `start` has ", length(st),
         " date(s) but `end` has ", length(en), ".", call. = FALSE)
  }
  if (any(is.na(st)) || any(is.na(en))) {
    stop("infant_vaccination.windows contains an unparseable date.",
         call. = FALSE)
  }
  if (any(en < st)) {
    i <- which(en < st)[1]
    stop("infant_vaccination.windows[", i, "] ends before it starts: ",
         format(st[i]), " to ", format(en[i]), call. = FALSE)
  }

  if (length(st) > 1) {
    ord <- order(st)
    s <- st[ord]; e <- en[ord]
    for (i in seq_len(length(s) - 1)) {
      if (s[i + 1] <= e[i]) {
        stop("infant_vaccination.windows ", i, " and ", i + 1, " overlap (",
             format(s[i]), "-", format(e[i]), " and ",
             format(s[i + 1]), "-", format(e[i + 1]), ").",
             "\n  Uptake is set per window, so a birth date inside two of ",
             "them has no single uptake and both coverage and doses would ",
             "silently take the first match.",
             call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}


# The catch-up programme's own coherence checks.
#
# Its failure modes are all silent. A cohort that ages past every declared
# band simply stops being protected part-way through the run; bounds that
# disagree with the data are never matched. Either way the submission
# looks plausible and merely shows too little effect, so both are caught
# before any modelling starts.
validate_catchup_config <- function(cfg, raw) {

  cu <- cfg$catchup
  if (!any(cfg$scenarios_df$catchup_coverage > 0)) return(invisible(TRUE))

  am <- cu$age_months
  if (length(am) != 2 || any(is.na(am)) || am[1] < 0 || am[2] <= am[1]) {
    stop("catchup_vaccination.age_at_campaign_months must be [min, max) in ",
         "months with 0 <= min < max; got: ",
         paste(am, collapse = ", "), call. = FALSE)
  }

  if (nrow(cu$age_bounds) == 0) {
    stop("catchup_vaccination.age_bounds is empty, so the cohort belongs to ",
         "no age band and the campaign could have no effect.", call. = FALSE)
  }

  if (any(!is.finite(cu$age_bounds$max_mo))) {
    open <- cu$age_bounds$age_group[!is.finite(cu$age_bounds$max_mo)]
    stop("catchup_vaccination.age_bounds needs a finite upper bound for ",
         "every band: an open-ended band has no birth window, so the ",
         "cohort cannot be tracked through it.",
         "\n  Open-ended: ", paste(open, collapse = ", "), call. = FALSE)
  }

  weeks <- sort(unique(parse_dates_strict(
    raw$admissions$target_end_date, cfg$input_date_formats$weekly_counts,
    "hospitaladmissions.csv$target_end_date")))
  weeks <- weeks[weeks >= ymd(cfg$data_start)]
  if (length(weeks) == 0) return(invisible(TRUE))

  # The declared bands must carry the cohort all the way to the end of the
  # horizon. If they run out first, the cohort ages into a band nobody
  # declared and its protection silently disappears mid-run.
  first_dose    <- min(do.call(c, lapply(cu$campaigns, function(cp) cp$start)))
  span_mo       <- as.numeric(max(weeks) - first_dose) / 30.4375
  oldest_at_end <- am[2] + span_mo
  covered_to    <- max(cu$age_bounds$max_mo)

  if (covered_to < oldest_at_end) {
    stop("catchup_vaccination.age_bounds only reaches ", round(covered_to),
         " months, but the cohort is up to ", ceiling(oldest_at_end),
         " months old by the last modelled week (", format(max(weeks)), ").",
         "\n  The cohort would age into an undeclared band and silently ",
         "lose its protection part-way through the run.",
         "\n  Add the older band(s) to catchup_vaccination.age_bounds.",
         call. = FALSE)
  }

  invisible(TRUE)
}


validate_global_config <- function(cfg, raw) {
  validate_countries(cfg, raw)
  validate_data_start(cfg, raw)
  validate_horizon_anchor(cfg, raw)
  validate_infant_windows(cfg)
  validate_campaign_windows(cfg, raw)
  validate_catchup_config(cfg, raw)
  validate_adult_config(cfg)
  validate_scenarios(cfg)
  validate_programme_disjoint(cfg)
  validate_population_coverage(cfg, raw)
  invisible(TRUE)
}
