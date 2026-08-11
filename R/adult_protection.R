# ============================================================ #
# ADULT PROGRAMME - calendar campaign + waning ensemble
#
# Entirely independent of the infant programme. Where the infant model
# derives eligibility from birth cohorts and waning from age band, the
# adult model derives coverage from a calendar campaign and waning from
# months since dose, using one whole VE curve per Monte-Carlo sample.
#
# Coverage is ONE-OFF and CUMULATIVE: the product is a single dose, so
# nobody is re-vaccinated. Each campaign recruits from the
# not-yet-vaccinated, and protection conferred earlier keeps decaying
# rather than resetting.
# ============================================================ #


# Expand campaigns into per-week coverage increments.
#
# Returns tibble(week, delta) where delta is the fraction of the
# eligible population newly vaccinated in that week. sum(delta) equals
# total_coverage (subject to campaign windows overlapping the model's
# weekly grid).
build_campaign_schedule <- function(campaigns, total_coverage, target_weeks) {

  empty <- tibble(week = as.Date(character()), delta = numeric())
  if (length(campaigns) == 0 || total_coverage <= 0) return(empty)

  map_dfr(seq_along(campaigns), function(i) {

    cp <- campaigns[[i]]
    wk <- target_weeks[target_weeks >= cp$start & target_weeks <= cp$end]

    if (length(wk) == 0) {
      warning("Adult campaign ", i, " (", format(cp$start), " to ",
              format(cp$end), ") does not overlap any modelled week; ",
              "it will contribute no coverage.", call. = FALSE)
      return(empty)
    }

    if (!identical(cp$profile, "uniform")) {
      stop("Unsupported campaign profile: '", cp$profile,
           "'. Only \"uniform\" is implemented.", call. = FALSE)
    }

    # Uniform accrual: the campaign's share of total coverage spread
    # evenly across the weeks its window covers.
    tibble(week  = wk,
           delta = total_coverage * cp$share / length(wk))
  })
}


# Coverage and mean residual VE per (age band, week, sample).
#
# This is a convolution of the campaign uptake curve with the VE curve:
#
#   coverage(t)     = sum over vaccination weeks w <= t of  delta(w)
#   protection(t)   = sum over w <= t of  delta(w) * VE(t - w)
#   residual_ve(t)  = protection(t) / coverage(t)
#
# protection(t) is coverage x mean-residual-VE, i.e. exactly the
# "cov * prot" term the scenario arithmetic needs. Because expected
# admissions are LINEAR in VE, the coverage-weighted mean is exact
# rather than an approximation.
#
# Parameters:
#   target_weeks        – Date vector, the model's weekly grid
#   waning_curves       – from load_waning_curves(): (sample, month, ve)
#   campaigns           – cfg$adult$campaigns
#   total_coverage      – cumulative coverage to reach across all campaigns
#   ve_beyond_curve     – "zero" or "hold_last", past the curve's last month
#   eligible_age_groups – bands the programme covers
#
# Returns tibble(age_group, target_end_date, sample, coverage, residual_ve).
build_adult_protection <- function(target_weeks,
                                   waning_curves,
                                   campaigns,
                                   total_coverage,
                                   ve_beyond_curve     = "zero",
                                   eligible_age_groups) {

  target_weeks <- sort(unique(as.Date(target_weeks)))
  samples      <- sort(unique(waning_curves$sample))

  zero_grid <- crossing(age_group       = eligible_age_groups,
                        target_end_date = target_weeks,
                        sample          = samples) %>%
    mutate(coverage = 0, residual_ve = 0)

  if (length(eligible_age_groups) == 0) return(zero_grid[0, ])

  sched <- build_campaign_schedule(campaigns, total_coverage, target_weeks)
  if (nrow(sched) == 0) return(zero_grid)

  max_month <- max(waning_curves$month)

  # Every (vaccination week, target week) pair where the dose precedes
  # the target week. Elapsed months use the mean month length; the VE
  # curve is monthly while the model grid is weekly, so some rounding is
  # unavoidable here.
  pairs <- crossing(sched, target_end_date = target_weeks) %>%
    filter(target_end_date >= week) %>%
    mutate(elapsed_mo = floor(as.numeric(target_end_date - week) / 30.4375),
           month_idx  = pmin(elapsed_mo, max_month))

  prot <- pairs %>%
    left_join(waning_curves, by = c("month_idx" = "month"),
              relationship = "many-to-many") %>%
    mutate(ve = if (identical(ve_beyond_curve, "zero")) {
                  ifelse(elapsed_mo > max_month, 0, ve)
                } else {
                  ve   # hold_last: month_idx was already clamped
                }) %>%
    group_by(target_end_date, sample) %>%
    summarise(coverage   = sum(delta),
              protection = sum(delta * ve),
              .groups    = "drop") %>%
    mutate(residual_ve = ifelse(coverage > 0, protection / coverage, 0)) %>%
    select(target_end_date, sample, coverage, residual_ve)

  # Weeks before the first dose carry no coverage at all
  full <- crossing(target_end_date = target_weeks, sample = samples) %>%
    left_join(prot, by = c("target_end_date", "sample")) %>%
    mutate(coverage    = coalesce(coverage, 0),
           residual_ve = coalesce(residual_ve, 0))

  # Coverage is assumed uniform across the eligible bands
  crossing(age_group = eligible_age_groups, full) %>%
    select(age_group, target_end_date, sample, coverage, residual_ve)
}
