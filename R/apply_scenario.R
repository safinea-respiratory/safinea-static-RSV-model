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


# Apply an INFANT vaccination scenario to baseline Monte-Carlo samples.
#
# Steps:
#   1. Compute the fraction of each age cohort born in a vaccination window.
#   2. Draw one vaccine effectiveness (VE) per sample from N(IE_mean, IE_sd).
#   3. Back-calculate the no-vaccination counterfactual from the observed data
#      using vacc_uptake_baseline (the uptake already reflected in the data).
#   4. Re-apply the scenario uptake to produce yes / no / total strata.
#
# Parameters:
#   df                   – baseline_df from simulate_weekly_age_fixed_margins
#   IE_mean / IE_sd      – vaccine effectiveness distribution parameters
#   vacc_uptake          – scenario uptake to apply (0–1)
#   vacc_start/vacc_end  – Date vectors (one element per season)
#   vacc_uptake_baseline – uptake already embedded in the observed baseline
#   waning_df            – data.frame(age_group, waning); waning ∈ [0,1]
#                          where 1 = full initial protection, 0 = none left.
#                          Bands absent from waning_df default to 0.
#   age_bounds           – data.frame(age_group, min_mo, max_mo)
apply_scenario <- function(df,
                           IE_mean, IE_sd,
                           vacc_uptake,
                           vacc_start, vacc_end,
                           vacc_uptake_baseline = 0,
                           waning_df,
                           age_bounds) {

  df %>%
    add_prop_born_in_window(vacc_start, vacc_end, age_bounds) %>%
    # One VE draw per sample (shared across all age groups in that draw)
    group_by(sample) %>%
    mutate(vacc_IE = rnorm(1, mean = IE_mean, sd = IE_sd)) %>%
    ungroup() %>%
    # Split each row into vaccinated (yes) and unvaccinated (no) proportions
    mutate(yes = vacc_uptake * prop_born_in_window,
           no  = 1 - yes) %>%
    rename(value_total = value) %>%
    pivot_longer(c("yes", "no"), values_to = "proportion", names_to = "immunisation") %>%
    left_join(waning_df, by = "age_group") %>%
    # A band with no waning entry is one this programme does not reach
    mutate(waning = coalesce(waning, 0)) %>%
    # Back-calculate the no-vaccination counterfactual from the observed data
    mutate(denom = 1 - (vacc_uptake_baseline * prop_born_in_window) *
                       (vacc_IE * waning)) %>%
    check_backcalc_denominator() %>%
    mutate(value_total_no_vax = value_total / denom) %>%
    select(-denom) %>%
    # Apply scenario: vaccinated stratum gets the (1 − VE × waning) reduction;
    # unvaccinated stratum remains at full risk (vacc_IE set to 0 for "no" rows)
    mutate(vacc_IE = ifelse(immunisation == "yes", vacc_IE, 0),
           value   = (1 - vacc_IE * waning) * proportion * value_total_no_vax) %>%
    select(-proportion, -vacc_IE) %>%
    pivot_wider(names_from = "immunisation", values_from = "value") %>%
    group_by(target_end_date, age_group, sample) %>%
    mutate(total = yes + no) %>%
    ungroup() %>%
    pivot_longer(c("yes", "no", "total"), values_to = "value", names_to = "immunisation")
}
