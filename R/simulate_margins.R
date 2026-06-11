# Sample weekly-by-age admission allocations with fixed margins.
#
# For each 4-week burden period independently, draws n_draws contingency
# tables of (weeks x ages) where:
#   - row sums match the observed weekly totals, and
#   - column sums match the observed per-age 4-week totals.
# Sampling is done via stats::r2dtable (Patefield algorithm).
#
# When the two marginal grand totals disagree (common when weekly and
# monthly data come from different sources or aggregation windows), the
# function reconciles by scaling one side before sampling.
#
# Depends on: round_preserve_sum() from R/utils.R
simulate_weekly_age_fixed_margins <- function(weekly_df, age_df,
                                              n_draws   = 200,
                                              seed      = 42,
                                              reconcile = c("scale_weekly",
                                                            "scale_annual")) {

  reconcile <- match.arg(reconcile)
  stopifnot(all(c("target_end_date", "weekly_rsv_hospitalisations") %in% names(weekly_df)))
  stopifnot(all(c("burden_start_date", "burden_end_date", "age_group",
                  "total_rsv_hospitalisations") %in% names(age_df)))
  set.seed(seed)

  wk <- weekly_df %>%
    select(target_end_date, weekly_total = weekly_rsv_hospitalisations) %>%
    arrange(target_end_date)

  bur <- age_df %>%
    select(burden_start_date, burden_end_date, age_group,
           age_total = total_rsv_hospitalisations) %>%
    arrange(burden_start_date, age_group)

  if (nrow(wk) == 0 || nrow(bur) == 0) return(tibble())

  period_starts <- bur %>% distinct(burden_start_date) %>% arrange(burden_start_date)

  map_dfr(seq_len(nrow(period_starts)), function(i) {

    this_start <- period_starts$burden_start_date[i]
    this_end   <- bur %>%
      filter(burden_start_date == this_start) %>%
      distinct(burden_end_date) %>%
      pull(burden_end_date)

    wk_ct <- wk %>%
      filter(target_end_date >= this_start, target_end_date <= this_end) %>%
      arrange(target_end_date)

    an_ct <- bur %>%
      filter(burden_start_date == this_start) %>%
      arrange(age_group)

    if (nrow(wk_ct) == 0 || nrow(an_ct) == 0) return(tibble())

    weeks <- wk_ct$target_end_date
    ages  <- an_ct$age_group

    # Reconcile grand totals when weekly and monthly sources disagree
    total_weekly <- sum(wk_ct$weekly_total, na.rm = TRUE)
    total_annual <- sum(an_ct$age_total,    na.rm = TRUE)
    if (total_weekly <= 0 || total_annual <= 0) return(tibble())

    if (abs(total_weekly - total_annual) > 1e-6) {
      if (reconcile == "scale_weekly") {
        wk_ct <- wk_ct %>% mutate(weekly_total = weekly_total * (total_annual / total_weekly))
      } else {
        an_ct <- an_ct %>% mutate(age_total = age_total * (total_weekly / total_annual))
      }
    }

    row_sums <- as.integer(round_preserve_sum(wk_ct$weekly_total))
    col_sums <- as.integer(round_preserve_sum(an_ct$age_total))

    # Independent rounding of each marginal can leave a ±1 unit discrepancy
    # in the grand totals; absorb it into the age column with the largest
    # fractional remainder.
    diff <- sum(row_sums) - sum(col_sums)
    if (diff != 0) {
      adj_idx <- order((an_ct$age_total - col_sums),
                       decreasing = (diff < 0))[1:abs(diff)]
      col_sums[adj_idx] <- col_sums[adj_idx] + sign(diff)
    }

    # r2dtable requires both dimensions >= 2; degenerate 1-row / 1-column
    # cases collapse to a deterministic vector replicated n_draws times.
    tabs <- if (length(col_sums) == 1) {
      replicate(n_draws, matrix(row_sums, ncol = 1), simplify = FALSE)
    } else if (length(row_sums) == 1) {
      replicate(n_draws, matrix(col_sums, nrow = 1), simplify = FALSE)
    } else {
      r2dtable(n = n_draws, r = row_sums, c = col_sums)
    }

    map_dfr(seq_len(n_draws), function(d) {
      M <- tabs[[d]]   # matrix [n_weeks x n_ages]
      tibble(
        target_end_date = ymd(rep(weeks, times = length(ages))),
        age_group       = rep(ages, each = length(weeks)),
        sample          = d,
        value           = as.vector(M)
      ) %>%
        group_by(target_end_date, sample) %>%
        mutate(weekly_total = sum(value)) %>%
        ungroup()
    })
  })
}
