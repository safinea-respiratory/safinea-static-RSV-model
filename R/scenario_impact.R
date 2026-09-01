# ============================================================ #
# Scenario impact summaries
#
# Reduces the sample-level submission to the standard impact tables:
# admissions averted, relative change, and doses, each expressed in
# absolute terms and per 100k of total or eligible population.
#
# Every interval is computed on PAIRED samples. Sample d of a scenario
# and sample d of the baseline share the same underlying Monte-Carlo
# draw of the (week x age) table, so the difference is taken within a
# draw and the quantiles are over those differences. Summarising each
# arm separately first and then subtracting would discard that pairing
# and inflate the interval.
#
# Output schema matches the ODE model's:
#   iso, scen, season, median, lo, hi, n_sim
# with lo/hi the 5th and 95th percentiles (a 90% interval).
# ============================================================ #


# Map every modelled week to the season whose burden window contains it.
#
# The season is DERIVED from the window dates rather than read from the
# `season` column, so this works on the upstream single-season files too.
# Where the column is present it is cross-checked, because a mis-shifted
# row in the generated two-season data would otherwise be invisible.
#
# Returns tibble(target_end_date, season).
week_season_map <- function(cfg, raw, weeks) {

  fmt   <- cfg$input_date_formats$burden_agegroups
  weeks <- sort(unique(as.Date(weeks)))

  win <- raw$burden %>%
    filter(country == cfg$countries_df$name[1]) %>%
    mutate(start = parse_dates_strict(start_date, fmt, "burden start_date"),
           end   = parse_dates_strict(end_date,   fmt, "burden end_date")) %>%
    distinct(start, end, .keep_all = TRUE) %>%
    mutate(derived = paste0(format(start, "%Y"), "/", format(end, "%Y"))) %>%
    arrange(start)

  if ("season" %in% names(win)) {
    bad <- win %>% filter(season != derived)
    if (nrow(bad) > 0) {
      stop("The burden file's `season` column disagrees with its window dates: ",
           paste0(bad$season, " vs ", bad$derived, collapse = "; "),
           call. = FALSE)
    }
  }

  # Overlapping windows would put a week in two seasons and count it twice.
  if (nrow(win) > 1 && any(win$start[-1] <= win$end[-nrow(win)])) {
    stop("Burden windows overlap, so a week could fall in two seasons.",
         call. = FALSE)
  }

  m <- crossing(target_end_date = weeks,
                win %>% select(start, end, season = derived)) %>%
    filter(target_end_date >= start, target_end_date <= end) %>%
    select(target_end_date, season)

  absent <- setdiff(as.character(weeks), as.character(m$target_end_date))
  if (length(absent) > 0) {
    stop(length(absent), " modelled week(s) fall outside every burden window, ",
         "so they belong to no season: ",
         paste(utils::head(absent, 5), collapse = ", "), call. = FALSE)
  }
  m
}


# The age bands a scenario actually vaccinates. A programme with zero
# uptake contributes nothing, so a scenario with no vaccination at all
# (the baseline) has an empty set and is excluded from the impact tables.
scenario_bands <- function(cfg, i) {
  sc <- cfg$scenarios_df
  c(if (sc$infant_uptake[i]  > 0) cfg$infant$age_bounds$age_group else character(0),
    if (sc$adult_coverage[i] > 0) sc$adult_age_groups[[i]]        else character(0))
}


# median + 90% interval of a paired quantity, per (location, scenario)
# and optionally per season.
#
# by_season = FALSE collapses the whole horizon into one figure and
# leaves the season column blank, which is how the ODE output reports
# doses.
summarise_paired <- function(df, value_col, n_sim, by_season = TRUE) {

  g <- if (by_season) c("location", "scenario_id", "season")
       else            c("location", "scenario_id")

  out <- df %>%
    group_by(across(all_of(g))) %>%
    summarise(median = median(.data[[value_col]], na.rm = TRUE),
              lo     = quantile(.data[[value_col]], 0.05, na.rm = TRUE),
              hi     = quantile(.data[[value_col]], 0.95, na.rm = TRUE),
              .groups = "drop")

  if (!by_season) out$season <- ""

  out %>%
    transmute(iso = location, scen = scenario_id, season,
              median, lo, hi, n_sim = n_sim) %>%
    arrange(scen, season, iso)
}


# All eight impact tables.
#
# Returns a named list; names become the CSV filenames in
# write_scenario_impact().
compute_scenario_impact <- function(submission, cfg, raw) {

  n_sim <- cfg$n_draws
  sc    <- cfg$scenarios_df
  wks   <- week_season_map(cfg, raw, unique(submission$target_end_date))

  # Scenarios that actually vaccinate someone. The baseline is the
  # reference everything is differenced against, so it is not itself a row.
  bands_by_scen <- setNames(lapply(seq_len(nrow(sc)), function(i) scenario_bands(cfg, i)),
                            sc$id)
  active   <- names(bands_by_scen)[lengths(bands_by_scen) > 0]
  baseline <- setdiff(sc$id, active)
  if (length(baseline) != 1) {
    stop("Impact tables need exactly one scenario with no vaccination to act ",
         "as the reference; found ", length(baseline),
         if (length(baseline)) paste0(" (", paste(baseline, collapse = ", "), ")"),
         ".", call. = FALSE)
  }
  union_bands <- sort(unique(unlist(bands_by_scen)))

  # ---- per-band admissions, per sample, tagged with season ----
  adm <- submission %>%
    filter(target == "rsv_hospitalisations",
           grepl("_immTotal$", pop_group), !grepl("^total_", pop_group)) %>%
    mutate(age_group = sub("_immTotal$", "", pop_group)) %>%
    select(location, scenario_id, age_group, output_type_id,
           target_end_date, value) %>%
    inner_join(wks, by = "target_end_date")

  # ---- populations ----
  pop <- raw$population %>%
    inner_join(cfg$countries_df %>% select(name, iso2), by = c("country" = "iso2")) %>%
    select(location = country, age_group, population)
  pop_total <- pop %>% group_by(location) %>%
    summarise(total_pop = sum(population), .groups = "drop")
  pop_for <- function(bands) {
    pop %>% filter(age_group %in% bands) %>%
      group_by(location) %>% summarise(elig_pop = sum(population), .groups = "drop")
  }

  # Paired scenario-vs-baseline admissions over a given age set, WITHIN
  # each season. Pairing is on (location, season, sample).
  paired_for <- function(band_fn) {
    map_dfr(active, function(sid) {
      b <- band_fn(sid)
      if (length(b) == 0) return(tibble())
      tot <- adm %>% filter(age_group %in% b) %>%
        group_by(location, scenario_id, season, output_type_id) %>%
        summarise(v = sum(value), .groups = "drop")
      inner_join(
        tot %>% filter(scenario_id == sid),
        tot %>% filter(scenario_id == baseline) %>%
          select(location, season, output_type_id, base = v),
        by = c("location", "season", "output_type_id")
      ) %>%
        mutate(averted = base - v,
               pct     = 100 * (v - base) / base)   # negative = averted
    })
  }

  by_scen  <- paired_for(function(sid) bands_by_scen[[sid]])
  by_union <- paired_for(function(sid) union_bands)

  # ---- doses, per sample ----
  # Dose rows are age-stratified and carry an "undefined" all-ages row
  # alongside the bands, so summing everything would double-count. Take
  # the all-ages row.
  #
  # Doses are NOT split by season: with a single campaign they all fall
  # in one season, and a second season of zeros would be noise rather
  # than information. This also matches the ODE output, which leaves the
  # season column blank for doses.
  doses <- submission %>%
    filter(target == "administered_doses", scenario_id %in% active,
           pop_group == "undefined") %>%
    group_by(location, scenario_id, output_type_id) %>%
    summarise(doses = sum(value), .groups = "drop")

  elig_scen  <- map_dfr(active, ~ pop_for(bands_by_scen[[.x]]) %>%
                          mutate(scenario_id = .x))
  elig_union <- pop_for(union_bands)

  add_pops <- function(df, elig) {
    df %>% left_join(pop_total, by = "location") %>%
      left_join(elig, by = intersect(names(elig), names(df)))
  }

  list(
    scenario_impact_1_pct_averted_scenario_ages =
      summarise_paired(by_scen,  "pct", n_sim),

    scenario_impact_2_pct_averted_union_ages =
      summarise_paired(by_union, "pct", n_sim),

    scenario_impact_3_abs_averted =
      summarise_paired(by_scen,  "averted", n_sim),

    scenario_impact_4a_doses_per100k_total =
      summarise_paired(add_pops(doses, elig_scen) %>%
                         mutate(x = doses / total_pop * 1e5), "x", n_sim,
                       by_season = FALSE),

    scenario_impact_4b_doses_per100k_eligible =
      summarise_paired(add_pops(doses, elig_scen) %>%
                         mutate(x = doses / elig_pop * 1e5), "x", n_sim,
                       by_season = FALSE),

    scenario_impact_5_averted_per100k_total =
      summarise_paired(add_pops(by_scen, elig_scen) %>%
                         mutate(x = averted / total_pop * 1e5), "x", n_sim),

    scenario_impact_6_averted_per100k_eligible =
      summarise_paired(add_pops(by_scen, elig_scen) %>%
                         mutate(x = averted / elig_pop * 1e5), "x", n_sim),

    scenario_impact_7_averted_per100k_union =
      summarise_paired(by_union %>% left_join(pop_total, by = "location") %>%
                         left_join(elig_union, by = "location") %>%
                         mutate(x = averted / elig_pop * 1e5), "x", n_sim)
  )
}


# Expected relative reduction, per scenario AND season, as a reference band.
#
# Over a scenario's own eligible bands the arithmetic is exactly
#   pct = -100 x coverage x residual_ve
# so the only question is which residual_ve to anchor the reference at.
# The band spans the range of dose ages that season actually contains:
#
#   lower edge (most reduction)   -100 x coverage x VE(youngest dose)
#   upper edge (least reduction)  -100 x coverage x VE(oldest dose)
#   dashed line                   the midpoint
#
# Each season spans ONE YEAR of waning: season 1 covers VE months 0-12,
# season 2 months 12-24, and so on by season index. So the second
# season's band sits well above the first's purely because its doses are
# a year older - which is the whole point of running two seasons off a
# single campaign.
#
# The bounds are deliberately the round 12-month marks rather than the
# exact range of dose ages the season contains. The exact range is
# slightly narrower (roughly 0-10 months for an autumn campaign), and
# using it would push points outside the band for a reason that has
# nothing to do with waning: admissions occurring BEFORE the campaign
# takes effect get no reduction at all and drag the season total down.
# The 12-month framing absorbs that.
#
# VE is the MEAN over ALL reps in the waning file, so the band tracks the
# central waning trajectory. Its width is WANING; the per-country error
# bars on the plot are Monte-Carlo uncertainty. Different quantities.
expected_reduction <- function(cfg, week_season) {

  curves <- read.csv(cfg$adult$waning_curves_path)
  ve_col <- curves[[cfg$adult$ve_target]]
  max_m  <- max(curves$month)

  mean_ve <- function(m) {
    m <- min(m, max_m)
    v <- ve_col[curves$month == m]
    if (length(v) == 0) {
      stop("The waning file has no month ", m, ", needed for the expected ",
           "reduction band.", call. = FALSE)
    }
    mean(v)
  }

  # Seasons ordered by when they start, so the index gives each one its
  # 12-month slice of the waning curve.
  ordered <- week_season %>%
    group_by(season) %>%
    summarise(first = min(target_end_date), .groups = "drop") %>%
    arrange(first) %>%
    mutate(idx = row_number() - 1L)

  sc <- cfg$scenarios_df

  map_dfr(seq_len(nrow(sc)), function(i) {
    cv <- sc$adult_coverage[i]
    if (cv <= 0) return(tibble())

    map_dfr(seq_len(nrow(ordered)), function(k) {
      lo_m <- 12L * ordered$idx[k]
      hi_m <- 12L * (ordered$idx[k] + 1L)
      tibble(scen = sc$id[i], season = ordered$season[k],
             lo_month = lo_m, hi_month = hi_m,
             lo       = -100 * cv * mean_ve(lo_m),   # freshest, most reduction
             hi       = -100 * cv * mean_ve(hi_m),   # oldest, least reduction
             expected = -100 * cv * (mean_ve(lo_m) + mean_ve(hi_m)) / 2)
    })
  })
}


# Write the impact tables as CSVs, one per list element.
write_scenario_impact <- function(impact, dir = "output/3_results") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- vapply(names(impact), function(nm) {
    p <- file.path(dir, paste0(nm, ".csv"))
    write.csv(impact[[nm]], p, row.names = FALSE, na = "")
    p
  }, character(1))
  message("Wrote ", length(paths), " impact tables to ", dir)
  invisible(paths)
}
