# ============================================================ #
# INFANT PROGRAMME - maternal / birth-dose product
#
# Everything in this file belongs to the infant programme only. Its
# defining property is that protection is conferred at (or before)
# birth, so time-since-vaccination equals age, and both eligibility and
# waning can be derived from the age band alone.
#
# The adult programme is entirely independent: eligibility comes from a
# calendar campaign and waning is indexed by months since dose, drawn
# from an ensemble. It shares no parameters with the code below.
# ============================================================ #


# For each (target_end_date, age_group) row, compute the fraction of
# the implied birth cohort that was born inside the union of vaccination
# windows. Used to split observed admissions between the "would have
# been vaccinated" and "would not" strata.
#
# age_bounds: data.frame(age_group, min_mo, max_mo) from the config,
#   giving each band's age span in months, [min inclusive, max exclusive).
#   Bands absent from age_bounds are ones the infant programme never
#   reached; they get prop_born_in_window = 0.
#
# vacc_start / vacc_end: Date vectors (one element per season)
add_prop_born_in_window <- function(df, vacc_start, vacc_end, age_bounds,
                                    age_col  = "age_group",
                                    date_col = "target_end_date") {

  bounds <- age_bounds %>% rename(!!age_col := age_group)

  # prop_born_in_window is a deterministic function of (age band, week)
  # alone - it does not depend on the Monte-Carlo sample. Computing it on
  # the distinct grid and joining back avoids repeating the rowwise work
  # once per draw, which is an n_draws-fold saving (100x here).
  key <- df %>% distinct(across(all_of(c(age_col, date_col))))

  key_prop <- key %>%
    left_join(bounds, by = age_col) %>%
    rowwise() %>%
    mutate(
      # Birth window: the range of birth dates for people aged
      # [min_mo, max_mo) months at the target week. NA bounds mean the
      # band is outside the infant programme entirely.
      birth_start = if (is.na(max_mo) || is.infinite(max_mo)) as.Date(NA)
                    else (!!sym(date_col)) %m-% months(max_mo),
      birth_end   = if (is.na(min_mo)) as.Date(NA)
                    else (!!sym(date_col)) %m-% months(min_mo),

      # Overlap between the cohort's birth window and each vaccination
      # season window, summed across seasons.
      inter_days = case_when(
        is.na(birth_start) | is.na(birth_end) ~ 0,
        birth_end <= birth_start              ~ 0,
        TRUE ~ {
          inter_starts <- pmax(birth_start, vacc_start)
          inter_ends   <- pmin(birth_end,   vacc_end)
          sum(pmax(0, as.numeric(inter_ends - inter_starts)), na.rm = TRUE)
        }
      ),
      denom_days = case_when(
        is.na(birth_start) | is.na(birth_end) ~ Inf,
        birth_end <= birth_start              ~ 0,
        TRUE                                  ~ as.numeric(birth_end - birth_start)
      ),
      prop_born_in_window = case_when(
        is.infinite(denom_days) | denom_days <= 0 ~ 0,
        TRUE ~ pmax(0, pmin(1, inter_days / denom_days))
      )
    ) %>%
    ungroup() %>%
    select(all_of(c(age_col, date_col)), prop_born_in_window)

  df %>% left_join(key_prop, by = c(age_col, date_col))
}


# Guard the no-vaccination back-calculation.
#
# The counterfactual divides observed admissions by
# (1 - uptake x prop_born x VE x waning). At the current parameters that
# bottoms out near 0.22, but an uptake and waning both close to 1 with a
# high VE draw would divide by ~0 and emit Inf admissions into the
# submission with no error anywhere.
# Two thresholds, because there are two distinct failure modes:
#   floor - at or below zero the counterfactual is mathematically
#           undefined and emits Inf or negative admissions. Hard error.
#   warn  - a small positive denominator is valid arithmetic but implies
#           an implausible inflation factor (1/denom). 0.05 corresponds
#           to scaling observed admissions up 20-fold.
check_backcalc_denominator <- function(df, floor = 1e-6, warn_below = 0.05) {

  worst <- suppressWarnings(min(df$denom, na.rm = TRUE))
  if (!is.finite(worst)) return(df)

  if (worst < floor) {
    stop("The no-vaccination back-calculation divides by ",
         "(1 - uptake x prop_born x VE x waning), which reached ",
         signif(worst, 3), ".",
         "\n  At or below zero the counterfactual is undefined and would ",
         "produce Inf or negative admissions.",
         "\n  Reduce infant_vaccination.baseline_uptake or the ",
         "waning_by_band values.",
         call. = FALSE)
  }

  if (worst < warn_below) {
    warning("The no-vaccination back-calculation divides by as little as ",
            signif(worst, 3), ", inflating observed admissions up to ",
            round(1 / worst), "-fold.",
            "\n  Check infant_vaccination.baseline_uptake and waning_by_band ",
            "- this counterfactual is unlikely to be plausible.",
            call. = FALSE)
  }

  df
}


# INFANT protection table.
#
# Emits the same two-column contract as build_adult_protection():
#   coverage    – fraction of the cell vaccinated
#   residual_ve – protection still held by a vaccinated individual
#
# For this product coverage comes from the birth-cohort overlap with the
# vaccination windows, and residual VE from the age band (age being the
# same thing as time-since-dose here).
#
# ve_draws: data.frame(sample, vacc_IE), drawn once per country and shared
#   across scenarios so a given sample index means the same VE world in
#   every scenario.
build_infant_protection <- function(grid, uptake, vacc_start, vacc_end,
                                    age_bounds, waning_df, ve_draws) {

  # Emit rows ONLY for bands this programme covers. Returning a row for
  # every band would collide with the adult table on bind_rows(), giving
  # two rows per key and silently multiplying admissions in the join.
  grid %>%
    filter(age_group %in% age_bounds$age_group) %>%
    add_prop_born_in_window(vacc_start, vacc_end, age_bounds) %>%
    left_join(waning_df, by = "age_group") %>%
    # A band with no waning entry is one this programme does not reach
    mutate(waning = coalesce(waning, 0)) %>%
    left_join(ve_draws, by = "sample") %>%
    mutate(coverage    = uptake * prop_born_in_window,
           residual_ve = vacc_IE * waning) %>%
    select(age_group, target_end_date, sample, coverage, residual_ve)
}


# A protection table must carry at most one row per
# (age_group, target_end_date, sample). Duplicates mean two programmes
# claimed the same band, and a left_join would multiply admissions
# instead of failing.
check_protection_unique <- function(protection, label) {
  key  <- c("age_group", "target_end_date", "sample")
  dups <- protection %>%
    count(across(all_of(key))) %>%
    filter(n > 1)

  if (nrow(dups) > 0) {
    bands <- sort(unique(dups$age_group))
    stop("The ", label, " protection table has duplicate rows for ",
         nrow(dups), " (age_group, week, sample) key(s).",
         "\n  Affected age group(s): ", paste(bands, collapse = ", "),
         "\n  Each age band must be claimed by at most one programme; a ",
         "duplicate would multiply admissions in the join.",
         call. = FALSE)
  }
  invisible(TRUE)
}


# Apply a vaccination scenario to baseline Monte-Carlo samples.
#
# Programme-agnostic: it consumes protection tables and knows nothing
# about birth cohorts, campaigns or waning curves. Each programme derives
# its own (coverage, residual_ve) independently; because the config
# forces the two to cover disjoint age bands, they can simply be
# bind_rows()-ed together before being passed in.
#
# Bands covered by no programme are absent from the tables and default to
# zero coverage, so they pass through at their observed values.
#
# The arithmetic:
#   no_vax = observed / (1 - coverage_base x residual_ve_base)
#   yes    = (1 - residual_ve) x coverage       x no_vax
#   no     =                     (1 - coverage) x no_vax
#   total  = yes + no
#
# Parameters:
#   df                  – baseline_df from simulate_weekly_age_fixed_margins
#   protection          – scenario protection table (age_group,
#                         target_end_date, sample, coverage, residual_ve)
#   protection_baseline – same shape, at the coverage already embedded in
#                         the observed data. Drives the back-calculation.
apply_scenario <- function(df, protection, protection_baseline) {

  key <- c("age_group", "target_end_date", "sample")

  # A duplicated key would silently multiply admissions in the joins
  # below rather than erroring - the two programmes must not both claim
  # a band. Cheap to check, and the failure is otherwise invisible.
  check_protection_unique(protection,          "scenario")
  check_protection_unique(protection_baseline, "baseline")

  df %>%
    rename(value_total = value) %>%
    left_join(protection, by = key) %>%
    left_join(protection_baseline %>%
                rename(coverage_base = coverage, residual_ve_base = residual_ve),
              by = key) %>%
    mutate(across(c(coverage, residual_ve, coverage_base, residual_ve_base),
                  ~ coalesce(.x, 0))) %>%
    # Back-calculate the no-vaccination counterfactual from the observed data
    mutate(denom = 1 - coverage_base * residual_ve_base) %>%
    check_backcalc_denominator() %>%
    mutate(value_total_no_vax = value_total / denom) %>%
    # Stratify: the vaccinated fraction carries the (1 - residual VE)
    # reduction, the unvaccinated fraction stays at full risk
    mutate(yes = (1 - residual_ve) * coverage       * value_total_no_vax,
           no  =                     (1 - coverage) * value_total_no_vax,
           total = yes + no) %>%
    select(-denom, -coverage, -residual_ve, -coverage_base, -residual_ve_base) %>%
    pivot_longer(c("yes", "no", "total"),
                 values_to = "value", names_to = "immunisation")
}
