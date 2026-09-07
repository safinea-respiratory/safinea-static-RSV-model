# ============================================================ #
# CAMPAIGN PROGRAMMES - calendar campaign + waning ensemble
#
# Shared machinery for every programme whose defining feature is that
# protection is conferred ON A DATE rather than at birth, so waning is
# indexed by MONTHS SINCE DOSE.
#
# Two programmes use it, and they differ in ONE respect only - how
# eligibility is defined:
#
#   adult    fixed_bands   the eligible group IS a set of age bands, and
#                          stays that set. A 67-year-old is in "65-69"
#                          all season, so "who we vaccinated" and "who is
#                          in the band" are the same people throughout.
#
#   catch-up birth_cohort  the eligible group is an age range ON THE
#                          CAMPAIGN DATE, i.e. a birth cohort. Infants
#                          age out of their band within weeks, so the
#                          cohort MOVES across bands while the band
#                          refills with children who were never eligible.
#
# That distinction is the whole reason this file is not simply the adult
# module: a static band list cannot express a moving cohort. Concretely,
# for an under-6-months campaign the "0-2mo" band is 100 % covered on
# campaign day and 0 % three months later, while "6-11mo" - a band you
# would never have listed - goes the other way.
#
# The two also draw their waning from different places, because they are
# different vaccines:
#
#   adult     an ensemble of whole VE curves, one per Monte-Carlo sample,
#             from data/vaccine/waning_curves.csv.
#   catch-up  the INFANT product's vacc_IE and waning_by_band from the
#             YAML - the same vaccine as the birth dose - re-indexed onto
#             a months-since-dose clock by build_infant_waning_curve().
#
# Both arrive here as tibble(sample, month, ve), so the convolution below
# is identical for either.
#
# Coverage is ONE-OFF and CUMULATIVE: the product is a single dose, so
# nobody is re-vaccinated. Each campaign recruits from the
# not-yet-vaccinated, and protection conferred earlier keeps decaying
# rather than resetting.
# ============================================================ #


# ---- Eligibility rules --------------------------------------------------
#
# Each returns a small spec consumed by eligibility_weight() below. They
# exist so that call sites read as intent ("this campaign targets a birth
# cohort") rather than as a mode string threaded through the arguments.

# Eligibility by age band, constant over time. The adult case.
eligibility_fixed_bands <- function(age_groups) {
  list(mode = "fixed_bands", age_groups = unique(as.character(age_groups)))
}


# Eligibility by age ON THE DAY OF VACCINATION, i.e. a birth cohort.
#
#   age_months  [min, max) age in months at the moment of the dose.
#               c(0, 6) means "everyone under six months old that day".
#   age_bounds  data.frame(age_group, min_mo, max_mo) for every band the
#               cohort can occupy over the modelled horizon - not just
#               the ones it occupies on campaign day. An under-6-months
#               cohort ends the season in "6-11mo" and the next one in
#               "1-4", so all four infant bands belong here. A band left
#               out simply receives no protection, which is silent, so
#               validate_catchup_config() checks the span instead.
eligibility_birth_cohort <- function(age_months, age_bounds) {
  stopifnot(length(age_months) == 2, age_months[1] < age_months[2])
  if (any(!is.finite(age_bounds$min_mo)) || any(!is.finite(age_bounds$max_mo))) {
    stop("birth-cohort eligibility needs finite month bounds for every ",
         "band it can reach; an open-ended band (max_mo = Inf) has no ",
         "birth window and cannot be tracked.", call. = FALSE)
  }
  list(mode       = "birth_cohort",
       lo_mo      = as.numeric(age_months[1]),
       hi_mo      = as.numeric(age_months[2]),
       age_bounds = age_bounds)
}


# Every band this rule can ever touch.
eligibility_bands <- function(el) {
  if (identical(el$mode, "fixed_bands")) el$age_groups
  else                                   el$age_bounds$age_group
}


# The eligibility weight w(band, t, w): of the people in `band` at target
# week `t`, what fraction was eligible for a dose given in week `w`?
#
#   fixed_bands   1 for every listed band, at every t. The band IS the
#                 eligible group, so the weight never changes and this
#                 reduces to the original adult behaviour exactly.
#
#   birth_cohort  the overlap of two birth-date windows, as a fraction of
#                 the band's own width:
#                   band at t    born in [t - max_mo, t - min_mo)
#                   cohort at w  born in [w - hi_mo,  w - lo_mo)
#                 Both windows are calendar-month arithmetic (%m-%), so
#                 month lengths are respected rather than approximated.
#
# `pairs` carries one row per (dose week, target week) with the campaign
# increment; this returns it expanded by age band with a w_elig column,
# dropping combinations the cohort never reaches.
eligibility_weight <- function(el, pairs) {

  if (identical(el$mode, "fixed_bands")) {
    return(crossing(pairs, age_group = el$age_groups) %>% mutate(w_elig = 1))
  }

  crossing(pairs, el$age_bounds) %>%
    mutate(
      band_born_from = target_end_date %m-% months(max_mo),
      band_born_to   = target_end_date %m-% months(min_mo),
      coh_born_from  = week            %m-% months(el$hi_mo),
      coh_born_to    = week            %m-% months(el$lo_mo),
      overlap_days   = pmax(0, as.numeric(pmin(band_born_to, coh_born_to) -
                                          pmax(band_born_from, coh_born_from))),
      band_days      = as.numeric(band_born_to - band_born_from),
      w_elig         = ifelse(band_days > 0, overlap_days / band_days, 0)
    ) %>%
    filter(w_elig > 0) %>%
    select(week, delta, target_end_date, elapsed_mo, month_idx,
           age_group, w_elig)
}


# ---- Campaign schedule --------------------------------------------------

# Expand campaigns into per-week coverage increments.
#
# Returns tibble(week, delta) where delta is the fraction of the
# eligible population newly vaccinated in that week. sum(delta) equals
# total_coverage (subject to campaign windows overlapping the model's
# weekly grid).
build_campaign_schedule <- function(campaigns, total_coverage, target_weeks,
                                    label = "Adult") {

  empty <- tibble(week = as.Date(character()), delta = numeric())
  if (length(campaigns) == 0 || total_coverage <= 0) return(empty)

  map_dfr(seq_along(campaigns), function(i) {

    cp <- campaigns[[i]]
    wk <- target_weeks[target_weeks >= cp$start & target_weeks <= cp$end]

    # A HARD ERROR, not a warning. A campaign that matches no modelled
    # week silently zeroes the whole scenario - doses, coverage and
    # effect - and the run still completes with a full submission of
    # zeros. That surfaces much later as something unrelated-looking,
    # so it has to stop here.
    if (length(wk) == 0) {
      near <- target_weeks[order(abs(as.numeric(target_weeks - cp$start)))][1:2]
      stop(label, " campaign ", i, " (", format(cp$start), " to ",
           format(cp$end), ") matches no modelled week, so it would ",
           "deliver zero coverage.",
           "\n  Campaign windows are matched against week-ending dates, ",
           "which are always ", weekdays(target_weeks[1]), "s.",
           "\n  ", format(cp$start), " is a ", weekdays(cp$start), ".",
           "\n  Nearest valid dates: ",
           paste(format(sort(near)), collapse = ", "),
           call. = FALSE)
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


# ---- Protection ---------------------------------------------------------

# Coverage and mean residual VE per (age band, week, sample).
#
# This is a convolution of the campaign uptake curve with the VE curve,
# weighted by eligibility:
#
#   coverage(b,t)    = sum over vaccination weeks w <= t of
#                        delta(w) * w_elig(b,t,w)
#   protection(b,t)  = sum over w <= t of
#                        delta(w) * w_elig(b,t,w) * VE(t - w)
#   residual_ve(b,t) = protection(b,t) / coverage(b,t)
#
# protection(b,t) is coverage x mean-residual-VE, i.e. exactly the
# "cov * prot" term the scenario arithmetic needs. Because expected
# admissions are LINEAR in VE, the coverage-weighted mean is exact
# rather than an approximation - and for the same reason the eligibility
# weight can be folded straight into the same sum.
#
# With fixed_bands eligibility w_elig is identically 1, so both sums
# collapse to the plain campaign convolution and the adult programme's
# numbers are unchanged.
#
# Parameters:
#   target_weeks    – Date vector, the model's weekly grid
#   waning_curves   – from load_waning_curves(): (sample, month, ve)
#   campaigns       – list of campaign windows
#   total_coverage  – cumulative coverage to reach across all campaigns
#   ve_beyond_curve – "zero" or "hold_last", past the curve's last month
#   eligibility     – from eligibility_fixed_bands() / _birth_cohort()
#   label           – programme name, used in error messages
#
# Returns tibble(age_group, target_end_date, sample, coverage, residual_ve).
build_campaign_protection <- function(target_weeks,
                                      waning_curves,
                                      campaigns,
                                      total_coverage,
                                      ve_beyond_curve = "zero",
                                      eligibility,
                                      label = "Adult") {

  target_weeks <- sort(unique(as.Date(target_weeks)))
  samples      <- sort(unique(waning_curves$sample))
  bands        <- eligibility_bands(eligibility)

  zero_grid <- crossing(age_group       = bands,
                        target_end_date = target_weeks,
                        sample          = samples) %>%
    mutate(coverage = 0, residual_ve = 0)

  if (length(bands) == 0) return(zero_grid[0, ])

  sched <- build_campaign_schedule(campaigns, total_coverage, target_weeks, label)
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

  # Who, in each band, is protected by a dose given in week w
  pairs <- eligibility_weight(eligibility, pairs)
  if (nrow(pairs) == 0) return(zero_grid)

  prot <- pairs %>%
    left_join(waning_curves, by = c("month_idx" = "month"),
              relationship = "many-to-many") %>%
    mutate(ve = if (identical(ve_beyond_curve, "zero")) {
                  ifelse(elapsed_mo > max_month, 0, ve)
                } else {
                  ve   # hold_last: month_idx was already clamped
                }) %>%
    group_by(age_group, target_end_date, sample) %>%
    summarise(coverage   = sum(delta * w_elig),
              protection = sum(delta * w_elig * ve),
              .groups    = "drop") %>%
    mutate(residual_ve = ifelse(coverage > 0, protection / coverage, 0)) %>%
    select(age_group, target_end_date, sample, coverage, residual_ve)

  # Bands and weeks the campaign has not reached carry no coverage at all
  zero_grid %>%
    select(age_group, target_end_date, sample) %>%
    left_join(prot, by = c("age_group", "target_end_date", "sample")) %>%
    mutate(coverage    = coalesce(coverage, 0),
           residual_ve = coalesce(residual_ve, 0))
}


# Backwards-compatible entry point for the adult programme.
#
# Kept because "adult protection" is what the diagnostics and the config
# talk about; it is a fixed-bands campaign and nothing more.
build_adult_protection <- function(target_weeks, waning_curves, campaigns,
                                   total_coverage, ve_beyond_curve = "zero",
                                   eligible_age_groups) {
  build_campaign_protection(
    target_weeks, waning_curves, campaigns, total_coverage, ve_beyond_curve,
    eligibility = eligibility_fixed_bands(eligible_age_groups),
    label       = "Adult")
}


# Re-read the INFANT waning table as a months-since-dose curve.
#
# infant_vaccination.waning_by_band is written per age band, but for a
# birth dose age and months-since-dose are the same number - a "6-11mo"
# infant is one 6 to 11 months past its dose. So that table already IS a
# months-since-dose curve; only its labelling says otherwise.
#
# The catch-up gives the SAME PRODUCT on a later date, which separates
# the two clocks: a child dosed at four months old sits in "6-11mo" two
# months later while being only two months post-dose. Re-indexing by
# month is what lets the shared campaign machinery apply the infant
# product's own waning on the right clock.
#
# Returns tibble(sample, month, ve) - the same shape load_waning_curves()
# produces for the adult ensemble, so build_campaign_protection() cannot
# tell the two apart. Here the between-sample spread comes from the
# vacc_IE draw rather than from an ensemble of whole curves.
build_infant_waning_curve <- function(age_bounds, waning_df, ve_draws) {

  finite <- age_bounds[is.finite(age_bounds$max_mo), ]
  if (nrow(finite) == 0) {
    stop("The infant age_bounds define no finite band, so no ",
         "months-since-dose curve can be built.", call. = FALSE)
  }

  months <- seq(0, max(finite$max_mo) - 1)

  # Which band's waning applies m months after the dose
  band_at <- vapply(months, function(m) {
    hit <- finite$age_group[finite$min_mo <= m & m < finite$max_mo]
    if (length(hit) == 0) NA_character_ else hit[1]
  }, character(1))

  tibble(month = months, age_group = band_at) %>%
    left_join(waning_df, by = "age_group") %>%
    # A month no band covers has no protection defined; treat as none
    mutate(waning = coalesce(waning, 0)) %>%
    crossing(ve_draws) %>%
    mutate(ve = vacc_IE * waning) %>%
    select(sample, month, ve) %>%
    arrange(sample, month)
}


# Entry point for the infant catch-up programme.
#
# Same machinery, birth-cohort eligibility, and the infant product's own
# waning re-indexed by months since dose. `age_months` is the age range
# eligible on the day of the dose; `age_bounds` covers every band the
# cohort can age into before the horizon ends.
#
# ve_beyond_curve is fixed at "zero": the curve runs to the end of the
# oldest infant band, past which the product confers nothing.
build_catchup_protection <- function(target_weeks, infant_waning, campaigns,
                                     total_coverage, age_months, age_bounds) {
  build_campaign_protection(
    target_weeks, infant_waning, campaigns, total_coverage,
    ve_beyond_curve = "zero",
    eligibility     = eligibility_birth_cohort(age_months, age_bounds),
    label           = "Catch-up")
}
