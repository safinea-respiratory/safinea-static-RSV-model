# ============================================================ #
# Submission output: validation and writing
#
# The input side is guarded by R/validate.R. This is the same
# fail-fast principle turned around to face the output: nothing else
# checks that what we are about to write is internally coherent, so a
# malformed submission would only be discovered downstream.
#
# Every check here is derived from the config and the data - no external
# schema is consulted.
# ============================================================ #


# Validate the assembled submission.
#
# Collects ALL failures and reports them together rather than stopping at
# the first, so one run tells you everything that is wrong.
#
# sample_draws limits the three PER-CELL ARITHMETIC invariants - dose
# totals, immYes + immNo == immTotal, and total_* == the sum over bands -
# to that many Monte-Carlo draws. Everything structural still runs on the
# whole table: columns, NAs, identifiers, horizon, pop_group membership,
# duplicate rows and grid completeness.
#
# Why that split is safe. The three sampled checks assert identities that
# format_submission() constructs row by row, so they hold for every draw
# or for none - a draw is not a place where they could fail selectively.
# They are also, together, almost the entire cost: they group by
# (location, scenario, horizon, draw, band), which at 28 countries x 16
# scenarios x 104 weeks x 100 draws is 51 million groups and takes over a
# quarter of an hour. Checking five draws is 20x less work.
#
# Set sample_draws = NULL (or >= n_draws) to check every draw.
validate_submission <- function(submission, cfg, sample_draws = 5) {

  problems <- character(0)
  note <- function(...) problems <<- c(problems, paste0(...))

  # NA-safe maximum. Without this an NA anywhere in `value` would make
  # the comparisons below evaluate to NA, and `if (NA)` is an R error -
  # the validator would crash instead of reporting the NA it had already
  # detected, along with everything else wrong.
  safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) 0 else max(x)
  }

  # ---- columns ----
  required <- c("round_id", "scenario_id", "target", "pop_group", "horizon",
                "target_end_date", "output_type", "output_type_id", "value",
                "location")
  absent <- setdiff(required, names(submission))
  if (length(absent) > 0) {
    # Nothing else can be trusted without the columns, so stop here.
    stop("Submission is missing required column(s): ",
         paste(absent, collapse = ", "),
         "\n  Columns present: ", paste(names(submission), collapse = ", "),
         call. = FALSE)
  }

  # ---- the subset the arithmetic invariants run on ----
  # Draws are picked evenly across the range rather than taking the first
  # few, so a fault confined to late draws is still in scope.
  all_ids <- unique(submission$output_type_id)
  n_all   <- length(all_ids)
  use_all <- is.null(sample_draws) || sample_draws >= n_all

  if (use_all) {
    arith    <- submission
    n_used   <- n_all
  } else {
    ord      <- all_ids[order(suppressWarnings(as.numeric(all_ids)), all_ids)]
    pick     <- ord[unique(round(seq(1, n_all, length.out = sample_draws)))]
    arith    <- submission %>% filter(output_type_id %in% pick)
    n_used   <- length(pick)
    message("validate_submission: per-cell arithmetic checked on ", n_used,
            " of ", n_all, " draws (sample_draws = ", sample_draws,
            "). Structural checks use every row.")
  }
  # Appended to those three findings so a count is never mistaken for a
  # count over the whole submission.
  scope <- if (use_all) "" else paste0(" [of ", n_used, " draw(s) checked]")

  # ---- no missing or nonsensical values ----
  for (col in required) {
    n_na <- sum(is.na(submission[[col]]))
    if (n_na > 0) note("column `", col, "` has ", n_na, " NA value(s)")
  }
  if (any(!is.finite(submission$value))) {
    note(sum(!is.finite(submission$value)), " non-finite value(s) in `value`")
  }
  neg <- sum(submission$value < -1e-9, na.rm = TRUE)
  if (neg > 0) note(neg, " negative value(s) in `value`")

  # ---- identifiers match the config ----
  if (!identical(sort(unique(submission$round_id)), sort(cfg$round_id))) {
    note("round_id is ",
         paste0('"', unique(submission$round_id), '"', collapse = ", "),
         "; config says \"", cfg$round_id, "\"")
  }
  chk_set <- function(col, expected, label) {
    got <- sort(unique(submission[[col]]))
    if (!setequal(got, expected)) {
      note(label, " mismatch.\n      in submission: ",
           paste(got, collapse = ", "),
           "\n      in config:     ", paste(sort(expected), collapse = ", "))
    }
  }
  chk_set("scenario_id", cfg$scenarios_df$id,  "scenario_id")
  chk_set("location",    cfg$countries_df$iso2, "location")

  if (!setequal(unique(submission$output_type), "sample")) {
    note("output_type must be \"sample\"; found ",
         paste(unique(submission$output_type), collapse = ", "))
  }
  if (!setequal(unique(submission$output_type_id), as.character(seq_len(cfg$n_draws)))) {
    note("output_type_id should be \"1\"..\"", cfg$n_draws, "\"; found ",
         length(unique(submission$output_type_id)), " distinct value(s)")
  }

  # ---- horizon is consistent with the dates ----
  bad_h <- submission %>%
    mutate(expected = as.integer((target_end_date - cfg$anchor) / 7)) %>%
    filter(horizon != expected)
  if (nrow(bad_h) > 0) {
    note(nrow(bad_h), " row(s) where horizon disagrees with ",
         "(target_end_date - anchor)/7")
  }

  # ---- pop_group belongs to the right target ----
  # Doses are age-stratified: each row is an age band, plus one
  # "undefined" all-ages row per (scenario, week). Hospitalisations use
  # <band>_<immStatus>, whose aggregate is "total_<immStatus>" - note the
  # two targets use different names for their aggregate row.
  # unique() BEFORE the regex: there are ~36 distinct pop_groups but tens
  # of millions of rows, and sub() over every row costs 4x what it needs to.
  data_bands <- submission %>%
    filter(target == "rsv_hospitalisations", grepl("_imm", pop_group)) %>%
    pull(pop_group) %>% unique() %>%
    sub("_imm(Yes|No|Total)$", "", .) %>%
    unique() %>% setdiff("total")

  d_bad <- submission %>%
    filter(target == "administered_doses",
           !pop_group %in% c(data_bands, "undefined"))
  if (nrow(d_bad) > 0) {
    note(nrow(d_bad), " administered_doses row(s) whose pop_group is ",
         "neither an age band nor \"undefined\": ",
         paste(head(unique(d_bad$pop_group), 5), collapse = ", "))
  }
  h_bad <- submission %>%
    filter(target == "rsv_hospitalisations", !grepl("_imm", pop_group))
  if (nrow(h_bad) > 0) {
    note(nrow(h_bad), " rsv_hospitalisations row(s) whose pop_group is not ",
         "<band>_<immStatus>: ",
         paste(head(unique(h_bad$pop_group), 5), collapse = ", "))
  }

  # ---- dose totals equal the sum over dose bands ----
  dose_chk <- arith %>%
    filter(target == "administered_doses") %>%
    mutate(is_total = pop_group == "undefined") %>%
    group_by(location, scenario_id, horizon, output_type_id) %>%
    summarise(d = abs(sum(value[is_total]) - sum(value[!is_total])),
              .groups = "drop")
  if (safe_max(dose_chk$d) > 1e-6) {
    note("administered_doses \"undefined\" rows do not equal the sum over age ",
         "bands in ", sum(dose_chk$d > 1e-6, na.rm = TRUE), " cell(s)", scope,
         "; worst ", signif(safe_max(dose_chk$d), 3))
  }

  # ---- no duplicated rows ----
  # select(all_of(...)) rather than [, key] because `submission` is a
  # data.table, where the latter does not mean what it means for a
  # data.frame.
  key <- c("location", "scenario_id", "target", "pop_group", "horizon",
           "output_type_id")
  n_dup <- nrow(submission) -
           nrow(distinct(submission %>% select(all_of(key))))
  if (n_dup > 0) {
    note(n_dup, " duplicate row(s) on (",
         paste(key, collapse = ", "), ")")
  }

  hosp <- arith %>% filter(target == "rsv_hospitalisations")

  # ---- immYes + immNo == immTotal ----
  strata <- hosp %>%
    filter(!grepl("^total_", pop_group)) %>%
    mutate(band = sub("_imm(Yes|No|Total)$", "", pop_group),
           st   = sub("^.*_imm", "", pop_group)) %>%
    group_by(location, scenario_id, horizon, output_type_id, band) %>%
    summarise(d = abs(sum(value[st == "Yes"]) + sum(value[st == "No"]) -
                        sum(value[st == "Total"])),
              .groups = "drop")
  if (safe_max(strata$d) > 1e-6) {
    note("immYes + immNo != immTotal in ",
         sum(strata$d > 1e-6, na.rm = TRUE),
         " cell(s)", scope, "; worst discrepancy ", signif(safe_max(strata$d), 3))
  }

  # ---- all-ages totals equal the sum over bands ----
  totals <- hosp %>%
    mutate(is_total = grepl("^total_", pop_group),
           st       = sub("^.*_imm", "", pop_group)) %>%
    group_by(location, scenario_id, horizon, output_type_id, st) %>%
    summarise(d = abs(sum(value[is_total]) - sum(value[!is_total])),
              .groups = "drop")
  if (safe_max(totals$d) > 1e-6) {
    note("total_* rows do not equal the sum over age bands in ",
         sum(totals$d > 1e-6, na.rm = TRUE), " cell(s)", scope, "; worst ",
         signif(safe_max(totals$d), 3))
  }

  # ---- the grid is complete ----
  for (tgt in unique(submission$target)) {
    sub_t <- submission %>% filter(target == tgt)
    expected <- length(unique(sub_t$location)) *
                length(unique(sub_t$scenario_id)) *
                length(unique(sub_t$horizon)) *
                length(unique(sub_t$output_type_id)) *
                length(unique(sub_t$pop_group))
    if (nrow(sub_t) != expected) {
      note("target \"", tgt, "\" has ", nrow(sub_t), " rows but a complete ",
           "grid would be ", expected,
           " (", nrow(sub_t) - expected, " difference)")
    }
  }

  if (length(problems) > 0) {
    stop("Submission failed validation (", length(problems), " problem(s)):\n",
         paste0("  - ", problems, collapse = "\n"), call. = FALSE)
  }

  message("Submission validated: ", format(nrow(submission), big.mark = ","),
          " rows, ", length(unique(submission$location)), " location(s), ",
          length(unique(submission$scenario_id)), " scenario(s).")
  invisible(TRUE)
}


# Write the submission to parquet, creating the output directory if needed.
#
# The default filename is derived from round_id so successive rounds do
# not overwrite one another.
# zstd rather than gzip. Measured on this submission's shape, gzip costs
# 7.6 s per million rows against zstd's 0.55 s - 25 minutes versus under
# two for a full round - for 7 % less file. Both read back in arrow and
# in nanoparquet.
SUBMISSION_COMPRESSION <- "zstd"

write_submission <- function(submission, cfg, path = NULL,
                             compression = SUBMISSION_COMPRESSION) {

  if (is.null(path)) {
    path <- file.path("output", paste0(cfg$round_id, "_staticModel.parquet"))
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  arrow::write_parquet(submission, path, compression = compression)

  message("Wrote ", format(nrow(submission), big.mark = ","), " rows to ", path,
          " (", round(file.size(path) / 1024^2, 2), " MB)")
  invisible(path)
}


# Concatenate per-country parquet parts into the single dense file the
# hub expects, WITHOUT ever holding the whole thing.
#
# The submission is ~14 GB in memory at 28 countries, 16 scenarios and
# 100 draws, so binding the parts and writing once needs more RAM than
# most machines have - and it is the reason a parallel run was killed
# with "error reading from connection". Arrow writes row group by row
# group, so peak memory is one country.
#
# Parts must share a schema; a mismatch would otherwise be written as a
# ragged file that only fails when something tries to read it.
write_submission_parts <- function(parts, cfg, path = NULL,
                                   compression = SUBMISSION_COMPRESSION) {

  if (length(parts) == 0) {
    stop("No submission parts to concatenate.", call. = FALSE)
  }

  # The parts must be exactly one per configured country.
  #
  # Callers are meant to pass the paths run_all_countries() returned. The
  # tempting alternative - globbing the parts directory - would silently
  # mix runs, and the schema check below cannot see that: parts from a
  # different run have the SAME schema. Checking against the config is
  # what makes the mistake an error rather than a wrong file.
  missing <- parts[!file.exists(parts)]
  if (length(missing) > 0) {
    stop("Submission part(s) not found:
  ",
         paste(missing, collapse = "
  "), call. = FALSE)
  }

  got  <- tools::file_path_sans_ext(basename(parts))
  want <- cfg$countries_df$iso2
  if (!setequal(got, want) || anyDuplicated(got) > 0) {
    stop("The parts given do not match the configured countries - one per ",
         "country is expected.",
         if (length(setdiff(want, got)))
           paste0("
  Missing: ", paste(setdiff(want, got), collapse = ", ")),
         if (length(setdiff(got, want)))
           paste0("
  Unexpected: ", paste(setdiff(got, want), collapse = ", ")),
         if (anyDuplicated(got) > 0)
           paste0("
  Duplicated: ",
                  paste(unique(got[duplicated(got)]), collapse = ", ")),
         "
  Pass the paths run_all_countries() returned in $parts rather ",
         "than listing the directory; parts left by an earlier run look ",
         "identical and would be merged in silently.",
         call. = FALSE)
  }
  if (is.null(path)) {
    path <- file.path("output", paste0(cfg$round_id, "_staticModel.parquet"))
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  sch  <- arrow::open_dataset(parts[1], format = "parquet")$schema
  sink <- arrow::FileOutputStream$create(path)
  on.exit(try(sink$close(), silent = TRUE), add = TRUE)

  w <- arrow::ParquetFileWriter$create(
    sch, sink,
    properties = arrow::ParquetWriterProperties$create(
      column_names = names(sch), compression = compression))

  n <- 0
  for (p in parts) {
    tb <- arrow::read_parquet(p, as_data_frame = FALSE)
    if (!identical(names(tb$schema), names(sch))) {
      stop("Submission part ", basename(p), " has different columns from ",
           basename(parts[1]), ".",
           "
  ", basename(parts[1]), ": ", paste(names(sch), collapse = ", "),
           "
  ", basename(p), ": ", paste(names(tb$schema), collapse = ", "),
           call. = FALSE)
    }
    w$WriteTable(tb, chunk_size = 1e6)
    n <- n + tb$num_rows
  }
  w$Close(); sink$close()

  message("Wrote ", format(n, big.mark = ","), " rows from ", length(parts),
          " part(s) to ", path, " (", round(file.size(path) / 1024^2, 2), " MB)")
  invisible(path)
}


# The one check that is NOT per country: every configured country must
# have produced a part.
#
# Each country validates its own rows as it is built, against a config
# narrowed to itself - so grid completeness, duplicates and the strata
# identities are all already covered. What a per-country check cannot
# see is a country that silently produced nothing at all.
validate_submission_set <- function(results, cfg) {

  got  <- vapply(results, function(r) r$location, character(1), USE.NAMES = FALSE)
  want <- cfg$countries_df$iso2

  absent <- setdiff(want, got)
  extra  <- setdiff(got,  want)
  if (length(absent) > 0 || length(extra) > 0) {
    stop("The submission parts do not match the configured countries.",
         if (length(absent)) paste0("
  Missing: ", paste(absent, collapse = ", ")),
         if (length(extra))  paste0("
  Unexpected: ", paste(extra, collapse = ", ")),
         call. = FALSE)
  }
  if (anyDuplicated(got) > 0) {
    stop("More than one part for: ",
         paste(unique(got[duplicated(got)]), collapse = ", "), call. = FALSE)
  }

  message("All ", length(want), " countries produced a part (",
          format(sum(vapply(results, function(r) r$n_rows, numeric(1))),
                 big.mark = ","), " rows).")
  invisible(TRUE)
}
