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


# Season label from the burden window, e.g. "2026/2027".
season_label <- function(cfg, raw) {
  d <- raw$burden %>%
    filter(country == cfg$countries_df$name[1]) %>%
    slice(1)
  s <- parse_dates_strict(d$start_date, cfg$input_date_formats$burden_agegroups,
                          "burden start_date")
  e <- parse_dates_strict(d$end_date, cfg$input_date_formats$burden_agegroups,
                          "burden end_date")
  paste0(format(s, "%Y"), "/", format(e, "%Y"))
}


# The age bands a scenario actually vaccinates. A programme with zero
# uptake contributes nothing, so a scenario with no vaccination at all
# (the baseline) has an empty set and is excluded from the impact tables.
scenario_bands <- function(cfg, i) {
  sc <- cfg$scenarios_df
  c(if (sc$infant_uptake[i]  > 0) cfg$infant$age_bounds$age_group else character(0),
    if (sc$adult_coverage[i] > 0) sc$adult_age_groups[[i]]        else character(0))
}


# median + 90% interval of a paired quantity, per (location, scenario).
summarise_paired <- function(df, value_col, season, n_sim) {
  df %>%
    group_by(location, scenario_id) %>%
    summarise(median = median(.data[[value_col]], na.rm = TRUE),
              lo     = quantile(.data[[value_col]], 0.05, na.rm = TRUE),
              hi     = quantile(.data[[value_col]], 0.95, na.rm = TRUE),
              .groups = "drop") %>%
    transmute(iso = location, scen = scenario_id, season = season,
              median, lo, hi, n_sim = n_sim) %>%
    arrange(scen, iso)
}


# All eight impact tables.
#
# Returns a named list; names become the CSV filenames in
# write_scenario_impact().
compute_scenario_impact <- function(submission, cfg, raw) {

  season <- season_label(cfg, raw)
  n_sim  <- cfg$n_draws
  sc     <- cfg$scenarios_df

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

  # ---- per-band admissions, per sample ----
  adm <- submission %>%
    filter(target == "rsv_hospitalisations",
           grepl("_immTotal$", pop_group), !grepl("^total_", pop_group)) %>%
    mutate(age_group = sub("_immTotal$", "", pop_group)) %>%
    select(location, scenario_id, age_group, output_type_id, value)

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

  # Paired scenario-vs-baseline admissions over a given age set.
  paired_for <- function(band_fn) {
    map_dfr(active, function(sid) {
      b <- band_fn(sid)
      if (length(b) == 0) return(tibble())
      tot <- adm %>% filter(age_group %in% b) %>%
        group_by(location, scenario_id, output_type_id) %>%
        summarise(v = sum(value), .groups = "drop")
      inner_join(
        tot %>% filter(scenario_id == sid),
        tot %>% filter(scenario_id == baseline) %>%
          select(location, output_type_id, base = v),
        by = c("location", "output_type_id")
      ) %>%
        mutate(averted = base - v,
               pct     = 100 * (v - base) / base)   # negative = averted
    })
  }

  by_scen  <- paired_for(function(sid) bands_by_scen[[sid]])
  by_union <- paired_for(function(sid) union_bands)

  # ---- doses, per sample ----
  # Dose rows are age-stratified and carry a "total" row alongside the
  # bands, so summing everything would double-count. Take the total row.
  doses <- submission %>%
    filter(target == "administered_doses", scenario_id %in% active,
           pop_group == "total") %>%
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
      summarise_paired(by_scen,  "pct", season, n_sim),

    scenario_impact_2_pct_averted_union_ages =
      summarise_paired(by_union, "pct", season, n_sim),

    scenario_impact_3_abs_averted =
      summarise_paired(by_scen,  "averted", season, n_sim),

    # Doses are not season-specific in the ODE output; the column is kept
    # for schema compatibility but left blank.
    scenario_impact_4a_doses_per100k_total =
      summarise_paired(add_pops(doses, elig_scen) %>%
                         mutate(x = doses / total_pop * 1e5), "x", "", n_sim),

    scenario_impact_4b_doses_per100k_eligible =
      summarise_paired(add_pops(doses, elig_scen) %>%
                         mutate(x = doses / elig_pop * 1e5), "x", "", n_sim),

    scenario_impact_5_averted_per100k_total =
      summarise_paired(add_pops(by_scen, elig_scen) %>%
                         mutate(x = averted / total_pop * 1e5), "x", season, n_sim),

    scenario_impact_6_averted_per100k_eligible =
      summarise_paired(add_pops(by_scen, elig_scen) %>%
                         mutate(x = averted / elig_pop * 1e5), "x", season, n_sim),

    scenario_impact_7_averted_per100k_union =
      summarise_paired(by_union %>% left_join(pop_total, by = "location") %>%
                         left_join(elig_union, by = "location") %>%
                         mutate(x = averted / elig_pop * 1e5), "x", season, n_sim)
  )
}


# Expected relative reduction per scenario, as a reference band.
#
# Over a scenario's own eligible bands the arithmetic is exactly
#   pct = -100 x coverage x residual_ve
# so the only question is which residual_ve to anchor the reference at.
# The band spans the WANING RANGE over the first year:
#
#   lower edge (most reduction)   -100 x coverage x VE(month 0)
#   upper edge (least reduction)  -100 x coverage x VE(month 12)
#   dashed line                   the midpoint of the two
#
# It therefore reads as "the effect if every dose were fresh" down to
# "the effect if every dose were a year old". A real campaign lands
# between the two, because by the time any given admission occurs its
# doses span a range of ages.
#
# VE is the MEAN over ALL reps in the waning file - all 500, not just the
# n_draws the model samples - so the band describes the central waning
# trajectory rather than any single realisation. Its width is therefore
# WANING, not Monte-Carlo uncertainty; the per-country intervals on the
# plot carry that.
#
# The band depends only on the config and the waning file, so a single
# one serves every country.
expected_reduction <- function(cfg) {

  curves <- read.csv(cfg$adult$waning_curves_path)
  ve_col <- curves[[cfg$adult$ve_target]]

  ve_at <- function(m) {
    v <- ve_col[curves$month == m]
    if (length(v) == 0) {
      stop("The waning file has no month ", m, ", which the expected ",
           "reduction band needs.",
           "
  months present: ", min(curves$month), " .. ", max(curves$month),
           "
  file: ", cfg$adult$waning_curves_path, call. = FALSE)
    }
    mean(v)
  }

  ve0  <- ve_at(0)
  ve12 <- ve_at(12)

  sc <- cfg$scenarios_df
  map_dfr(seq_len(nrow(sc)), function(i) {
    cv <- sc$adult_coverage[i]
    if (cv <= 0) return(tibble())
    tibble(scen     = sc$id[i],
           lo       = -100 * cv * ve0,               # fresh dose
           hi       = -100 * cv * ve12,              # dose 12 months old
           expected = -100 * cv * (ve0 + ve12) / 2)  # midpoint
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
