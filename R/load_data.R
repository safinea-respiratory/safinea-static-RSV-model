# Parse a character date column with an explicit format, strictly.
#
# R's as.Date() is dangerous here: given "01/09/2025" with no format it
# returns 0001-09-20 rather than NA, so a mis-formatted input silently
# produces dates ~2000 years out and every downstream join quietly
# yields zero rows. This wrapper requires the format up front, rejects
# anything that fails to parse, and additionally rejects implausible
# years so that class of error can never pass silently again.
parse_dates_strict <- function(x, format, what, min_year = 1900, max_year = 2100) {

  d <- as.Date(as.character(x), format = format)

  if (any(is.na(d))) {
    bad <- unique(as.character(x)[is.na(d)])
    stop("Could not parse ", what, " with format '", format, "'.",
         "\n  Unparseable value(s): ",
         paste0('"', utils::head(bad, 5), '"', collapse = ", "),
         if (length(bad) > 5) " ..." else "",
         "\n  Fix the data, or correct input_date_formats in the config.",
         call. = FALSE)
  }

  yr <- as.integer(format(d, "%Y"))
  if (any(yr < min_year | yr > max_year)) {
    bad_i <- which(yr < min_year | yr > max_year)
    stop("Parsing ", what, " with format '", format, "' produced ",
         "implausible date(s) - the format is almost certainly wrong.",
         "\n  e.g. \"", as.character(x)[bad_i[1]], "\" -> ", format(d[bad_i[1]]),
         "\n  Expected years between ", min_year, " and ", max_year, ".",
         call. = FALSE)
  }

  d
}


# Rename source age-group labels using the optional alias map.
# Labels with no alias pass through unchanged.
apply_age_aliases <- function(x, aliases) {
  if (is.null(aliases) || length(aliases) == 0) return(x)
  m   <- unlist(aliases)
  out <- unname(m[x])
  ifelse(is.na(out), x, out)
}


# Load and parse the YAML configuration file.
#
# Returns the raw cfg list augmented with pre-parsed dates and two
# normalised sub-lists, cfg$infant and cfg$adult, so that callers never
# reach into the nested YAML structure or re-parse dates.
#
# The infant and adult programmes are kept strictly separate: they share
# no parameters, and each sub-list is self-contained.
load_config <- function(path = "config/static_model.yaml") {
  cfg <- yaml::read_yaml(path)

  cfg$anchor <- as.Date(cfg$submission_horizon_anchor)

  # Countries: full name for target-data, ISO2 for auxiliary-data and
  # for the submission's `location` column.
  cfg$countries_df <- data.frame(
    name = vapply(cfg$countries, function(x) as.character(x$name), character(1)),
    iso2 = vapply(cfg$countries, function(x) as.character(x$iso2), character(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )

  # ---- Infant programme (maternal / birth-dose) ----
  inf <- cfg$infant_vaccination

  # Waning indexed by AGE BAND, because for this product age equals
  # time since vaccination.
  waning_df <- data.frame(
    age_group = names(inf$waning_by_band),
    waning    = unlist(inf$waning_by_band, use.names = FALSE),
    row.names = NULL
  )

  # Month bounds per band, [min inclusive, max exclusive).
  age_bounds <- data.frame(
    age_group = names(inf$age_bounds),
    min_mo    = vapply(inf$age_bounds, function(b) as.numeric(b[[1]]), numeric(1)),
    max_mo    = vapply(inf$age_bounds, function(b) as.numeric(b[[2]]), numeric(1)),
    row.names = NULL
  )

  cfg$infant <- list(
    vacc_IE_mean    = inf$vacc_IE$mean,
    vacc_IE_sd      = inf$vacc_IE$sd,
    waning_df       = waning_df,
    age_bounds      = age_bounds,
    vacc_start      = ymd(unlist(inf$windows$start)),
    vacc_end        = ymd(unlist(inf$windows$end)),
    baseline_uptake = inf$baseline_uptake,

    # Which band infant doses are counted against. This is a birth-dose
    # product, so they are attributed to the youngest band - the birth
    # cohort itself. Note that for a purely MATERNAL product the dose
    # arguably belongs to the mother's band instead; the choice only
    # matters once infant uptake is non-zero.
    dose_age_group  = age_bounds$age_group[which.min(age_bounds$min_mo)]
  )

  # ---- Adult programme (calendar campaign) ----
  # No VE level parameter here by design: the waning ensemble carries
  # the uncertainty, one whole curve per Monte-Carlo sample.
  adu <- cfg$adult_vaccination

  # Campaigns: normalise dates and fill in equal shares when omitted.
  campaigns <- lapply(adu$campaigns, function(cp) {
    list(start   = as.Date(cp$start),
         end     = as.Date(cp$end),
         share   = if (is.null(cp$share)) NA_real_ else as.numeric(cp$share),
         profile = if (is.null(cp$profile)) "uniform" else cp$profile)
  })
  if (length(campaigns) > 0 && all(vapply(campaigns, function(c) is.na(c$share), logical(1)))) {
    eq <- 1 / length(campaigns)
    campaigns <- lapply(campaigns, function(c) { c$share <- eq; c })
  }

  cfg$adult <- list(
    eligible_age_groups = unlist(adu$eligible_age_groups),
    waning_curves_path  = adu$waning_curves,
    ve_target           = adu$ve_target,
    ve_beyond_curve     = adu$ve_beyond_curve,
    campaigns           = campaigns,
    baseline_coverage   = adu$baseline_coverage
  )

  # ---- Scenarios ----
  # One row per submitted scenario. Both programmes are set independently,
  # so any combination of infant and adult uptake is expressible.
  #
  # adult_age_groups is a list-column: a scenario may target its own set
  # of bands, falling back to the programme-wide default. This makes
  # targeting strategies comparable at fixed coverage (e.g. 65+ vs 75+).
  # A tibble is used rather than data.frame because it carries
  # list-columns cleanly.
  cfg$scenarios_df <- tibble::tibble(
    id             = vapply(cfg$scenarios, function(s) as.character(s$id), character(1)),
    infant_uptake  = vapply(cfg$scenarios, function(s) as.numeric(s$infant_uptake), numeric(1)),
    adult_coverage = vapply(cfg$scenarios, function(s) as.numeric(s$adult_coverage), numeric(1)),
    adult_age_groups = lapply(cfg$scenarios, function(s) {
      if (is.null(s$adult_age_groups)) cfg$adult$eligible_age_groups
      else unlist(s$adult_age_groups)
    })
  )

  cfg
}


# Every adult band referenced anywhere: the programme-wide default plus
# any per-scenario override. Validation and the dose denominator both
# need the union, not just the default.
all_adult_bands <- function(cfg) {
  sort(unique(c(cfg$adult$eligible_age_groups,
                unlist(cfg$scenarios_df$adult_age_groups))))
}


# Read every input CSV once.
#
# All four files are country-agnostic; the per-country loaders below just
# filter these frames. Reading them inside the per-country loop instead
# would mean 60+ full-file reads for a 28-country run.
#
# Returns a list of raw data frames, to be passed to the loaders and to
# validate_countries().
load_raw_inputs <- function(
    cfg,
    weekly_file = "data/epidemiological/hospitaladmissions.csv",
    burden_file = "data/epidemiological/hospitalburden_agegroups.csv") {

  read_required <- function(path, label) {
    if (is.null(path) || !file.exists(path)) {
      stop(label, " not found: ", if (is.null(path)) "<unset>" else path,
           call. = FALSE)
    }
    read.csv(path)
  }

  list(
    admissions = read_required(weekly_file, "Weekly admissions file"),
    burden     = read_required(burden_file, "Age burden file"),
    births     = read_required(cfg$births_file, "Births file"),
    population = read_required(cfg$population_file, "Population file"),
    labels     = list(admissions = basename(weekly_file),
                      burden     = basename(burden_file),
                      births     = basename(cfg$births_file),
                      population = basename(cfg$population_file))
  )
}


# Select one country's weekly admissions and age-stratified burden from
# the raw target-data frames.
#
# Counts are used directly: the burden file gives absolute admissions per
# band rather than proportions, so no reconstruction is needed.
#
# Returns a named list:
#   $admissions — weekly aggregate counts
#                 (target_end_date, weekly_rsv_hospitalisations)
#   $burden     — per-age counts over the burden window
#                 (burden_start_date, burden_end_date, age_group,
#                  total_rsv_hospitalisations)
#
# NOTE ON TIME RESOLUTION: the burden file gives ONE total per age band
# for the whole season, so the fixed-margin sampler draws a single
# (weeks x ages) table per season rather than one per 4-week period.
# The age mix is therefore constrained only across the season as a
# whole, not within it. That is an honest reflection of what this data
# constrains, but it produces wider per-week age uncertainty than
# period-level age data would.
#
# Age-group labels come from the data itself; the only transformation
# applied is the optional alias map from the config.
load_epidemiological_data <- function(cfg, country_name, raw) {

  fmt_weekly <- cfg$input_date_formats$weekly_counts
  fmt_burden <- cfg$input_date_formats$burden_agegroups

  raw_adm <- raw$admissions
  raw_bur <- raw$burden

  check_country_present(country_name, raw_adm$country, raw$labels$admissions)
  check_country_present(country_name, raw_bur$country, raw$labels$burden)

  admissions <- raw_adm %>%
    filter(country == country_name) %>%
    mutate(target_end_date = parse_dates_strict(
                               target_end_date, fmt_weekly,
                               "hospitaladmissions.csv$target_end_date")) %>%
    select(target_end_date, weekly_rsv_hospitalisations) %>%
    arrange(target_end_date) %>%
    setDT()

  burden <- raw_bur %>%
    filter(country == country_name) %>%
    mutate(burden_start_date = parse_dates_strict(
                                 start_date, fmt_burden,
                                 "hospitalburden_agegroups.csv$start_date"),
           burden_end_date   = parse_dates_strict(
                                 end_date, fmt_burden,
                                 "hospitalburden_agegroups.csv$end_date"),
           age_group         = apply_age_aliases(age_group,
                                                 cfg$age_group_aliases)) %>%
    select(burden_start_date, burden_end_date, age_group,
           total_rsv_hospitalisations) %>%
    setDT()

  list(admissions = admissions, burden = burden)
}


# Error with the available options listed if the configured country is
# absent from an input file.
check_country_present <- function(country, column, file_label) {
  available <- sort(unique(column))
  if (!country %in% available) {
    stop("Country \"", country, "\" not found in ", file_label, ".",
         "\n  Available: ", paste(available, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}


# Load monthly births for one country.
#
# RespiCompass auxiliary-data supplies births already mapped onto the
# 2026-09 .. 2027-08 scenario period for every country, so no
# forward-projection of a historical pattern is needed. Keyed by ISO2.
load_births_data <- function(cfg, country_iso2, raw) {

  check_country_present(country_iso2, raw$births$country, raw$labels$births)

  raw$births %>%
    filter(country == country_iso2) %>%
    mutate(date = parse_dates_strict(date, cfg$input_date_formats$births,
                                     paste0(raw$labels$births, "$date"))) %>%
    select(country, date, births) %>%
    arrange(date)
}


# Load population by country x age band. Keyed by ISO2.
#
# The denominator for the adult administered-doses track: weekly adult
# doses are that week's coverage increment times the combined population
# of the eligible bands.
load_population_data <- function(cfg, country_iso2, raw) {

  check_country_present(country_iso2, raw$population$country,
                        raw$labels$population)

  raw$population %>%
    filter(country == country_iso2) %>%
    mutate(age_group = apply_age_aliases(age_group, cfg$age_group_aliases)) %>%
    select(country, age_group, population)
}


# Load the ADULT vaccine waning ensemble.
#
# The source file holds one complete VE-over-time curve per `rep`. Each
# Monte-Carlo sample takes one entire curve, so within-sample
# correlation across months is preserved by construction and no
# distributional assumption is required. This is deliberately unlike the
# infant programme, whose uncertainty is a parametric draw on vacc_IE.
#
# Reps are independent realisations, so the first n_draws are used.
# Errors if n_draws exceeds the ensemble size.
#
# Returns a data.frame(sample, month, ve) where `sample` is 1..n_draws
# and `ve` is the column named by cfg$adult$ve_target.
load_waning_curves <- function(cfg, n_draws) {

  path      <- cfg$adult$waning_curves_path
  ve_target <- cfg$adult$ve_target

  if (!file.exists(path)) {
    stop("Adult waning curve file not found: ", path, call. = FALSE)
  }

  curves <- read.csv(path)

  required <- c("rep", "month", ve_target)
  absent   <- setdiff(required, names(curves))
  if (length(absent) > 0) {
    stop("Waning curve file is missing required column(s): ",
         paste(absent, collapse = ", "),
         "\n  file: ", path,
         "\n  columns found: ", paste(names(curves), collapse = ", "),
         call. = FALSE)
  }

  reps   <- sort(unique(curves$rep))
  n_reps <- length(reps)

  if (n_draws > n_reps) {
    stop("n_draws (", n_draws, ") exceeds the number of reps in the adult ",
         "waning ensemble (", n_reps, ").",
         "\n  Reduce n_draws, or supply an ensemble with more reps.",
         "\n  file: ", path,
         call. = FALSE)
  }

  # Every rep must supply a complete, gap-free curve starting at month 0
  per_rep <- curves %>%
    group_by(rep) %>%
    summarise(n_months   = n(),
              min_month  = min(month),
              n_distinct = n_distinct(month),
              .groups    = "drop")

  expected_months <- max(per_rep$n_months)
  bad <- per_rep %>%
    filter(min_month != 0 | n_months != n_distinct | n_months != expected_months)

  if (nrow(bad) > 0) {
    stop("Waning curves are malformed for rep(s): ",
         paste(utils::head(bad$rep, 10), collapse = ", "),
         if (nrow(bad) > 10) " ..." else "",
         "\n  Each rep must start at month 0 and carry the same number of ",
         "unique, non-duplicated months.",
         "\n  file: ", path,
         call. = FALSE)
  }

  ve_vals <- curves[[ve_target]]
  if (any(is.na(ve_vals)) || any(ve_vals < 0 | ve_vals > 1)) {
    stop("Column '", ve_target, "' contains NA or values outside [0, 1].",
         "\n  file: ", path, call. = FALSE)
  }

  keep <- reps[seq_len(n_draws)]

  curves %>%
    filter(rep %in% keep) %>%
    mutate(sample = match(rep, keep),
           ve     = .data[[ve_target]]) %>%
    select(sample, month, ve) %>%
    arrange(sample, month)
}
