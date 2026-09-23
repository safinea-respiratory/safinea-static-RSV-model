# Run the full pipeline for a single country.
#
# Everything downstream of data loading is country-agnostic, so this is
# a thin wrapper: load that country's data, sample, apply scenarios,
# format, and stamp the `location` column.
#
# Parameters:
#   cfg          – parsed config
#   country_name – full name, used to select rows from target-data
#   country_iso2 – two-letter code, used for auxiliary-data and emitted
#                  as the submission's `location`
#   adult_waning – the waning ensemble, loaded once and shared
#   quiet        – suppress the informational validation messages
#                  (they are identical for every country)
#
# Returns a list:
#   $submission      – hospitalisations + doses, with `location`
#   $baseline_df     – the raw Monte-Carlo samples
#   $adult_protection – named list of coverage/residual-VE per scenario
run_country <- function(cfg, country_name, country_iso2,
                        adult_waning, raw, quiet = FALSE) {

  epi    <- load_epidemiological_data(cfg, country_name, raw)
  validate_age_groups(cfg, epi, quiet = quiet)
  births <- load_births_data(cfg, country_iso2, raw)
  popn   <- load_population_data(cfg, country_iso2, raw)

  # ---- Baseline Monte-Carlo ----
  baseline_df <- simulate_weekly_age_fixed_margins(
    weekly_df = epi$admissions %>% filter(target_end_date   >= ymd(cfg$data_start)),
    age_df    = epi$burden     %>% filter(burden_start_date >= ymd(cfg$data_start)),
    n_draws   = cfg$n_draws,
    seed      = cfg$mc_seed
  )

  if (nrow(baseline_df) == 0) {
    n_wk  <- sum(epi$admissions$target_end_date   >= ymd(cfg$data_start))
    n_bur <- sum(epi$burden$burden_start_date >= ymd(cfg$data_start))
    stop("No baseline samples produced for ", country_name, ".",
         "\n  After filtering on data_start (", cfg$data_start, "): ",
         n_wk, " weekly row(s), ", n_bur, " burden row(s).",
         "\n  Unfiltered the data spans:",
         "\n    burden window:     ", format(min(epi$burden$burden_start_date)),
         " .. ", format(max(epi$burden$burden_end_date)),
         "\n    weekly admissions: ", format(min(epi$admissions$target_end_date)),
         " .. ", format(max(epi$admissions$target_end_date)),
         "\n  The sampler needs both sides non-empty.",
         call. = FALSE)
  }

  # ---- Protection tables ----
  # Both programmes emit the same (coverage, residual_ve) contract, and
  # the config forces them onto disjoint age bands, so the two tables can
  # simply be stacked. apply_scenario() then needs to know nothing about
  # birth cohorts or campaigns.
  grid  <- baseline_df %>% distinct(age_group, target_end_date, sample)
  weeks <- unique(baseline_df$target_end_date)

  # One VE draw per sample, shared across scenarios so a given sample
  # index means the same VE world in every arm.
  ve_draws <- tibble(
    sample  = sort(unique(baseline_df$sample)),
    vacc_IE = rnorm(length(unique(baseline_df$sample)),
                    mean = cfg$infant$vacc_IE_mean,
                    sd   = cfg$infant$vacc_IE_sd)
  )

  # The infant table's expensive part - the rowwise birth-window overlap -
  # does not depend on uptake, so it is built ONCE here and each scenario
  # just scales it. With 16 scenarios this is the difference between one
  # rowwise pass per country and seventeen.
  infant_base <- build_infant_base(grid,
                                   cfg$infant$vacc_start, cfg$infant$vacc_end,
                                   cfg$infant$age_bounds, cfg$infant$waning_df,
                                   ve_draws)

  # ---- Scenarios ----
  # Every scenario is an elementwise multiply of the baseline draw, so
  # this is done as [week, band, draw] array arithmetic rather than by
  # joining and pivoting a long frame once per scenario. Same rows, same
  # values; see R/scenario_engine.R for why it is exact.
  submission_pre <- scenario_submission(
    baseline_df  = baseline_df,
    cfg          = cfg,
    adult_waning = adult_waning,
    ve_draws     = ve_draws,
    infant_base  = infant_base
  )

  sc <- cfg$scenarios_df

  # ---- Administered doses ----
  doses_df <- build_dose_table(
    baseline_df     = baseline_df,
    births_df       = births,
    population_df   = popn,
    vacc_start      = cfg$infant$vacc_start,
    vacc_end        = cfg$infant$vacc_end,
    cfg             = cfg,
    round_id        = cfg$round_id,
    anchor          = cfg$anchor,
    output_type_ids = unique(submission_pre$output_type_id)
  )

  submission <- bind_rows(submission_pre, doses_df) %>%
    mutate(location = country_iso2)

  # Adult coverage / residual VE per scenario, kept for inspection via
  # plot_adult_protection(). These now feed the admissions above rather
  # than sitting unused.
  adult_protection <- setNames(
    lapply(seq_len(nrow(sc)), function(i) {
      build_adult_protection(weeks, adult_waning, cfg$adult$campaigns,
                             sc$adult_coverage[i], cfg$adult$ve_beyond_curve,
                             sc$adult_age_groups[[i]])
    }),
    sc$id
  )

  list(
    submission       = submission,
    baseline_df      = baseline_df,
    adult_protection = adult_protection
  )
}


# Run every country listed in the config and bind the results.
#
# Countries are independent: each gets its own data, its own sampler
# draws and its own submission rows, distinguished by `location`. The
# same mc_seed is used throughout, so the vaccine-effectiveness draws
# are shared across countries - defensible, since it is the same vaccine.
#
# workers > 1 runs them on a PSOCK cluster. Independence is what makes
# that safe: no country reads another's state, and each seeds its own
# sampler from cfg$mc_seed, so the result does not depend on how the work
# was divided or on how many workers ran it.
#
# keep_by_country = FALSE drops the per-country objects once their
# submission rows have been taken. They exist for the diagnostic plots
# and are worth several hundred MB per country at submission scale, so
# holding all 28 is usually the largest thing in the session.
# Free physical memory in GB, or NA where it cannot be read.
free_ram_gb <- function() {
  out <- tryCatch(
    switch(Sys.info()[["sysname"]],
      Windows = {
        # CIM rather than wmic: wmic is deprecated and absent on
        # current Windows 11, where it returns nothing at all.
        x <- system2("powershell",
                     c("-NoProfile", "-Command",
                       "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"),
                     stdout = TRUE, stderr = FALSE)
        as.numeric(trimws(x[length(x)])) / 1024^2     # KB -> GB
      },
      Linux = {
        x <- readLines("/proc/meminfo")
        as.numeric(gsub("[^0-9]", "",
                        grep("^MemAvailable", x, value = TRUE))) / 1024^2
      },
      NA_real_),
    error = function(e) NA_real_, warning = function(w) NA_real_)
  if (length(out) != 1 || !is.finite(out) || out <= 0) NA_real_ else out
}


# Rough size of one country's submission, in GB. Rows are the product of
# the grid; 76 bytes/row is what the assembled table measures at.
estimate_country_gb <- function(cfg, bytes_per_row = 76) {
  n_pg   <- length(unlist(cfg$age_group_order)) * 3 + 3      # bands + total_*
  n_dose <- length(all_adult_bands(cfg)) + 1                 # bands + undefined
  weeks  <- 104                                              # order of magnitude
  rows   <- nrow(cfg$scenarios_df) * weeks * cfg$n_draws * (n_pg + n_dose)
  rows * bytes_per_row / 1024^3
}


run_all_countries <- function(cfg, adult_waning, raw,
                              workers         = 1L,
                              keep_by_country = TRUE) {

  cs <- cfg$countries_df
  n  <- nrow(cs)

  one <- function(i) {
    r <- run_country(cfg, cs$name[i], cs$iso2[i], adult_waning, raw,
                     quiet = (i > 1))
    if (keep_by_country) r else list(submission = r$submission)
  }

  workers <- max(1L, min(as.integer(workers), n, parallel::detectCores()))

  # Memory, not cores, is what limits this. Each worker holds a whole
  # country's submission while it builds it - roughly 500 MB at 16
  # scenarios and 100 draws - and the parent is meanwhile accumulating
  # every country it has already been handed. Ask for more workers than
  # RAM allows and a worker is killed mid-flight, which surfaces as
  # "error reading from connection" and says nothing about memory.
  if (workers > 1L) {
    per_gb  <- estimate_country_gb(cfg)
    free_gb <- free_ram_gb()
    if (is.finite(free_gb)) {
      # The parent ends up holding every country; workers add their own
      # copy on top. Leave a third of free memory as headroom for the
      # copies rbindlist makes.
      budget <- free_gb * 0.66 - per_gb * n
      fit    <- max(1L, floor(budget / per_gb))
      if (fit < workers) {
        message("Limiting to ", fit, " worker(s): ", signif(free_gb, 3),
                " GB free, ~", signif(per_gb, 2), " GB per country, and the ",
                "result itself is ~", signif(per_gb * n, 3), " GB.")
        workers <- fit
      }
    }
  }

  if (workers == 1L) {
    results <- lapply(seq_len(n), function(i) {
      message("Running ", cs$name[i], " (", cs$iso2[i], ") [", i, "/", n, "]")
      one(i)
    })
  } else {
    message("Running ", n, " countries on ", workers, " workers")
    cl <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # PSOCK workers start empty: they need the package namespaces and the
    # model's own functions before they can run a country.
    parallel::clusterEvalQ(cl, {
      suppressMessages({
        library(dplyr); library(tidyr); library(purrr); library(tibble)
        library(lubridate); library(data.table); library(yaml)
      })
      for (f in c("utils", "validate", "simulate_margins", "apply_scenario",
                  "adult_protection", "scenario_engine", "load_data",
                  "format_submission", "build_doses", "run_model")) {
        source(file.path("R", paste0(f, ".R")))
      }
      NULL
    })
    parallel::clusterExport(cl, c("cfg", "cs", "adult_waning", "raw",
                                  "one", "keep_by_country"),
                            envir = environment())
    results <- parallel::parLapplyLB(cl, seq_len(n), one)
  }
  names(results) <- cs$iso2

  submission <- rbindlist(lapply(results, `[[`, "submission"))
  if (!keep_by_country) results <- NULL

  list(submission = submission, by_country = results)
}
