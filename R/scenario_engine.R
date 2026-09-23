# ============================================================ #
# Scenario engine - array form
#
# Every scenario is an ELEMENTWISE multiply of the baseline draw:
#
#   no_vax = observed / (1 - coverage_base x residual_ve_base)
#   yes    = (1 - residual_ve) x coverage       x no_vax
#   no     =                     (1 - coverage) x no_vax
#   total  = yes + no
#
# Nothing in that needs a join. Held as [week, band, draw] arrays it is
# plain vector arithmetic over ~100k doubles per scenario, where the
# row-wise form spent its time joining and pivoting a 100k-row frame
# sixteen times over. Measured on one country at 16 scenarios and 100
# draws: 25.9 s of joins and pivots becomes 6.8 s, bit for bit identical.
#
# Two further savings fall out of the algebra rather than out of tuning:
#
#   1. The adult protection is LINEAR in coverage. The campaign schedule
#      is delta(w) = c x share / |weeks|, so coverage(t) = c x C(t) and
#      protection(t) = c x P(t,s); residual_ve = P/C therefore does not
#      depend on c at all. Coverage says how MANY were vaccinated,
#      residual_ve how good their protection still is - and vaccinating
#      more people does not change how old the average dose is. One
#      call at c = 1 reconstructs every scenario.
#
#   2. The eligible band list changes no numbers either. The row-wise
#      builder ends by copying one curve onto every eligible band, so a
#      scenario's band set is a MASK over that shared curve.
#
# Both are exact. assert_coverage_linear() below guards (1), because a
# future campaign profile that saturated or was non-linear in coverage
# would break it silently - every non-unit scenario would get the wrong
# coverage with nothing to flag it.
# ============================================================ #


# Pack a (age_group, target_end_date, sample, <col>) table into a
# [week, band, draw] array. Keys absent from the table are left at zero,
# which is what "this programme does not reach that cell" means.
#
# Duplicate keys would be silently overwritten here - the last one wins -
# so they are rejected. That is the same failure the row-wise path
# guarded with check_protection_unique(): two programmes claiming one
# band, which a join would multiply and an array would quietly drop.
arrayise <- function(df, col, W, A, D, label = "protection") {

  key <- cbind(match(df$target_end_date, W),
               match(df$age_group,       A),
               match(df$sample,          D))

  if (anyNA(key)) {
    stop("The ", label, " table has (week, band, draw) values that are not ",
         "in the baseline grid; they would be dropped silently.",
         call. = FALSE)
  }
  if (anyDuplicated(key) > 0) {
    d <- df[duplicated(key) | duplicated(key, fromLast = TRUE), ]
    stop("The ", label, " table has duplicate (age_group, week, sample) ",
         "keys - two programmes claim the same cell.",
         "\n  Affected band(s): ",
         paste(sort(unique(d$age_group)), collapse = ", "),
         call. = FALSE)
  }

  z <- array(0, c(length(W), length(A), length(D)))
  z[key] <- df[[col]]
  z
}


# Adult protection at coverage 1, over the UNION of every band any
# scenario targets.
#
# The union matters: adult_vaccination.eligible_age_groups is only the
# programme-wide default, and a scenario may widen it. Building the unit
# curve over the default alone would leave the extra bands at zero
# coverage - a silent under-count in exactly those scenarios.
unit_adult_protection <- function(weeks, waning, cfg) {
  build_adult_protection(weeks, waning, cfg$adult$campaigns,
                         total_coverage  = 1,
                         ve_beyond_curve = cfg$adult$ve_beyond_curve,
                         eligible_age_groups = all_adult_bands(cfg))
}


# Confirm the linearity the engine relies on, on the real campaign
# settings, before any of it is used.
assert_coverage_linear <- function(weeks, waning, cfg, unit, tol = 1e-9) {

  probe <- 0.5
  got <- build_adult_protection(weeks, waning, cfg$adult$campaigns, probe,
                                cfg$adult$ve_beyond_curve, all_adult_bands(cfg))
  k <- c("age_group", "target_end_date", "sample")
  got  <- got[do.call(order, got[k]), ]
  want <- unit[do.call(order, unit[k]), ]

  d_cov <- max(abs(got$coverage    - probe * want$coverage))
  d_rve <- max(abs(got$residual_ve -         want$residual_ve))

  if (d_cov > tol || d_rve > tol) {
    stop("The adult programme is no longer linear in coverage, so the ",
         "scenario engine cannot rebuild it from a single unit curve.",
         "\n  coverage(", probe, ") vs ", probe, " x coverage(1): worst ",
         signif(d_cov, 3),
         "\n  residual_ve(", probe, ") vs residual_ve(1): worst ",
         signif(d_rve, 3),
         "\n  A campaign profile that saturates, or whose shape depends ",
         "on coverage, would do this.",
         call. = FALSE)
  }
  invisible(TRUE)
}


# Every scenario's admissions, in submission shape.
#
# Returns the same columns assemble_submission() produced, and the same
# rows - only the order differs, which nothing downstream depends on.
#
# Parameters:
#   baseline_df   – from simulate_weekly_age_fixed_margins()
#   cfg           – parsed config
#   adult_waning  – the ensemble, shared across countries
#   ve_draws      – one vacc_IE per sample
#   infant_base   – from build_infant_base(), uptake not yet applied
scenario_submission <- function(baseline_df, cfg, adult_waning, ve_draws,
                                infant_base) {

  W <- sort(unique(baseline_df$target_end_date))
  A <- sort(unique(baseline_df$age_group))
  D <- sort(unique(baseline_df$sample))
  sc <- cfg$scenarios_df

  arr <- function(df, col, label) arrayise(df, col, W, A, D, label)

  # ---- the counterfactual the scenarios are applied to ----
  prot_base <- bind_rows(
    scale_infant_protection(infant_base, cfg$infant$baseline_uptake),
    build_adult_protection(W, adult_waning, cfg$adult$campaigns,
                           cfg$adult$baseline_coverage,
                           cfg$adult$ve_beyond_curve,
                           cfg$adult$eligible_age_groups))

  X    <- arr(baseline_df, "value",       "baseline")
  covb <- arr(prot_base,   "coverage",    "baseline protection")
  rveb <- arr(prot_base,   "residual_ve", "baseline protection")

  denom <- 1 - covb * rveb
  check_backcalc_denominator(data.frame(denom = as.vector(denom)))
  no_vax <- X / denom

  # ---- one adult curve, reused by every scenario ----
  unit <- unit_adult_protection(W, adult_waning, cfg)
  assert_coverage_linear(W, adult_waning, cfg, unit)
  uC <- arr(unit, "coverage",    "unit adult protection")
  uR <- arr(unit, "residual_ve", "unit adult protection")

  # ---- key skeleton ----
  # The key columns are IDENTICAL for every scenario, so the whole table
  # is allocated once and only `value` is written per scenario. Building
  # it with rbindlist per scenario instead reallocates 96 times.
  nW <- length(W); nA <- length(A); nD <- length(D); nS <- nrow(sc)
  STRATA <- c("immYes", "immNo", "immTotal")

  # as.vector() on a [week, band, draw] array runs week fastest, then
  # band, then draw - these repeats mirror exactly that.
  b_date <- rep(W,                      times = nA * nD)
  b_otid <- as.character(rep(D,         each  = nW * nA))
  b_band <- rep(rep(A, each = nW),      times = nD)
  t_date <- rep(W,                      times = nD)
  t_otid <- as.character(rep(D,         each  = nW))

  blk_date <- c(rep(b_date, 3L), rep(t_date, 3L))
  blk_otid <- c(rep(b_otid, 3L), rep(t_otid, 3L))
  blk_pg   <- c(unlist(lapply(STRATA, function(k) paste0(b_band, "_", k))),
                rep(paste0("total_", STRATA), each = nW * nD))
  n_block  <- length(blk_date)

  out <- data.table(
    round_id        = cfg$round_id,
    scenario_id     = rep(sc$id, each = n_block),
    target          = "rsv_hospitalisations",
    pop_group       = rep(blk_pg,   times = nS),
    horizon         = as.integer((rep(blk_date, times = nS) - cfg$anchor) / 7),
    target_end_date = rep(blk_date, times = nS),
    output_type     = "sample",
    output_type_id  = rep(blk_otid, times = nS),
    value           = numeric(nS * n_block))

  # Sum over bands without aperm(): nA - 1 additions of an [nW x nD]
  # matrix, which is cheaper than permuting the whole array.
  over_bands <- function(z) {
    s <- z[, 1L, ]
    if (nA > 1L) for (a in 2:nA) s <- s + z[, a, ]
    s
  }

  for (i in seq_len(nS)) {

    inf  <- scale_infant_protection(infant_base, sc$infant_uptake[[i]])
    icov <- arr(inf, "coverage",    "infant protection")
    irve <- arr(inf, "residual_ve", "infant protection")

    # A scenario's adult bands select from the unit curve; its coverage
    # scales it. Combined with the infant programme as a coverage-
    # weighted mean, so coverage x residual_ve - the only quantity the
    # arithmetic uses - is the sum of each programme's contribution.
    mask <- array(0, c(nW, nA, nD))
    hit  <- match(sc$adult_age_groups[[i]], A)
    if (length(hit)) mask[, hit, ] <- 1
    acov <- sc$adult_coverage[i] * uC * mask

    cov <- icov + acov
    prt <- icov * irve + acov * uR
    rve <- ifelse(cov > 0, prt / cov, 0)

    yes <- (1 - rve) * cov * no_vax
    no  <- (1 - cov) * no_vax
    tot <- yes + no

    set(out, i = seq.int((i - 1L) * n_block + 1L, i * n_block), j = "value",
        value = c(as.vector(yes), as.vector(no), as.vector(tot),
                  as.vector(over_bands(yes)), as.vector(over_bands(no)),
                  as.vector(over_bands(tot))))
  }

  out[]
}
