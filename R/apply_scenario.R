# For each (target_end_date, age_group) row, compute the fraction of
# the implied birth cohort that was born inside the union of vaccination
# windows. Used to split observed admissions between the "would have
# been vaccinated" and "would not" strata.
#
# vacc_start / vacc_end: Date vectors (one element per season)
add_prop_born_in_window <- function(df, vacc_start, vacc_end,
                                    age_col  = "age_group",
                                    date_col = "target_end_date") {

  # Map age-band labels to [min_months, max_months) of cohort age at
  # target_end_date. 65+y has no finite upper bound.
  age_map <- tibble::tibble(
    !!age_col := c("0-2mo", "3-5mo", "6-11mo", "1-4y", "5-64y", "65+y"),
    min_mo   = c(0,  3,  6, 12,  60,  780),
    max_mo   = c(3,  6, 12, 60, 780,  Inf)
  )

  df %>%
    left_join(age_map, by = age_col) %>%
    rowwise() %>%
    mutate(
      # Birth window: the range of birth dates for children aged
      # [min_mo, max_mo) months at the target week.
      birth_start = if (is.infinite(max_mo)) as.Date(NA)
                    else (!!sym(date_col)) %m-% months(max_mo),
      birth_end   = (!!sym(date_col)) %m-% months(min_mo),

      # Overlap between the cohort's birth window and each vaccination
      # season window, summed across seasons.
      inter_days = case_when(
        is.na(birth_start)       ~ 0,
        birth_end <= birth_start ~ 0,
        TRUE ~ {
          inter_starts <- pmax(birth_start, vacc_start)
          inter_ends   <- pmin(birth_end,   vacc_end)
          sum(pmax(0, as.numeric(inter_ends - inter_starts)), na.rm = TRUE)
        }
      ),
      denom_days = case_when(
        is.na(birth_start)       ~ Inf,
        birth_end <= birth_start ~ 0,
        TRUE                     ~ as.numeric(birth_end - birth_start)
      ),
      prop_born_in_window = case_when(
        is.infinite(denom_days) | denom_days <= 0 ~ 0,
        TRUE ~ pmax(0, pmin(1, inter_days / denom_days))
      )
    ) %>%
    ungroup() %>%
    select(-min_mo, -max_mo, -birth_start, -birth_end, -inter_days, -denom_days)
}


# Apply a vaccination scenario to baseline Monte-Carlo samples.
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
#                          where 1 = no protection remaining, 0 = full
apply_scenario <- function(df,
                           IE_mean, IE_sd,
                           vacc_uptake,
                           vacc_start, vacc_end,
                           vacc_uptake_baseline = 0,
                           waning_df) {

  df %>%
    add_prop_born_in_window(vacc_start, vacc_end) %>%
    # One VE draw per sample (shared across all age groups in that draw)
    group_by(sample) %>%
    mutate(vacc_IE = rnorm(1, mean = IE_mean, sd = IE_sd)) %>%
    ungroup() %>%
    # Split each row into vaccinated (yes) and unvaccinated (no) proportions
    mutate(vacc_uptake = vacc_uptake,
           yes = vacc_uptake * prop_born_in_window,
           no  = 1 - yes) %>%
    rename(value_total = value) %>%
    pivot_longer(c("yes", "no"), values_to = "proportion", names_to = "immunisation") %>%
    left_join(waning_df, by = "age_group") %>%
    # Back-calculate the no-vaccination counterfactual from the observed data
    mutate(value_total_no_vax =
             value_total / (1 - (vacc_uptake_baseline * prop_born_in_window) *
                                (vacc_IE * waning))) %>%
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
