# Diagnostic plots for the static RSV model output.
#
# All functions return a ggplot object so callers can further customise
# (add themes, save with ggsave, etc.) before displaying.
#
# Call from launch.R after the submission tables are assembled, e.g.:
#   print(plot_baseline_samples(baseline_df))
#   ggsave("output/scenario_comparison.png", plot_scenario_comparison(submission_pre))


# ---- Plot 1: Baseline Monte-Carlo uncertainty -------------------------
# Shows the distribution of weekly total RSV admissions (all age groups
# summed) across all MC draws. Ribbons cover the 50% and 95% credible
# intervals; the line is the draw-wise median.
# Purpose: verify that the fixed-margin sampler produces a plausible
# spread and that the baseline data look reasonable before applying
# any vaccination correction.
plot_baseline_samples <- function(baseline_df) {
  baseline_df %>%
    group_by(target_end_date, sample) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    group_by(target_end_date) %>%
    summarise(
      median = median(total),
      lo95   = quantile(total, 0.025),
      hi95   = quantile(total, 0.975),
      lo50   = quantile(total, 0.25),
      hi50   = quantile(total, 0.75),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = target_end_date)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "steelblue", alpha = 0.20) +
    geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "steelblue", alpha = 0.35) +
    geom_line(aes(y = median), colour = "steelblue", linewidth = 0.8) +
    labs(
      title    = "Baseline: weekly RSV hospital admissions (all ages)",
      subtitle = "Shaded bands: 50 % and 95 % Monte-Carlo intervals",
      x = "Week ending", y = "Admissions"
    ) +
    theme_bw()
}


# ---- Plot 2: Scenario comparison ----------------------------------------
# Compares median total-immunisation weekly admissions across all scenarios
# with 95 % MC intervals. Useful for judging the scale of the vaccination
# effect (high_vacc vs. no_vacc) relative to the observed baseline.
# Uses the pop_group = "total_immTotal" rows, which sum across all age bands.
plot_scenario_comparison <- function(submission_df) {
  submission_df %>%
    filter(target    == "rsv_hospitalisations",
           pop_group == "total_immTotal") %>%
    group_by(scenario_id, target_end_date) %>%
    summarise(
      median = median(value),
      lo95   = quantile(value, 0.025),
      hi95   = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = target_end_date, colour = scenario_id, fill = scenario_id)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.8) +
    labs(
      title    = "Scenario comparison: total weekly RSV admissions",
      subtitle = "Median + 95 % MC interval, all-ages total (immTotal)",
      x = "Week ending", y = "Admissions",
      colour = "Scenario", fill = "Scenario"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}


# ---- Plot 3: Age-group breakdown for a single scenario ------------------
# Shows median weekly admissions disaggregated by age band for one
# scenario. Use this to check that age-disaggregation produces a
# clinically plausible profile (e.g., infant groups peak earlier and
# higher relative to their population size).
# scenario: one of "baseline", "no_vacc", "high_vacc"
plot_age_breakdown <- function(submission_df, scenario = "baseline") {
  age_order <- c("0-2mo", "3-5mo", "6-11mo", "1-4y", "5-64y", "65+y")

  # pop_group format is "<age>_immTotal" for age-specific total rows
  submission_df %>%
    filter(target      == "rsv_hospitalisations",
           scenario_id == scenario,
           grepl("_immTotal$", pop_group),
           !grepl("^total_", pop_group)) %>%
    mutate(age_group = sub("_immTotal$", "", pop_group)) %>%
    filter(age_group %in% age_order) %>%
    mutate(age_group = factor(age_group, levels = age_order)) %>%
    group_by(age_group, target_end_date) %>%
    summarise(median = median(value), .groups = "drop") %>%
    ggplot(aes(x = target_end_date, y = median, colour = age_group)) +
    geom_line(linewidth = 0.7) +
    labs(
      title    = paste0("Age breakdown — scenario: '", scenario, "'"),
      subtitle = "Median weekly RSV admissions by age group (immTotal)",
      x = "Week ending", y = "Admissions", colour = "Age group"
    ) +
    theme_bw() +
    theme(legend.position = "right")
}


# ---- Plot 4: Weekly administered dose schedule --------------------------
# Shows the weekly number of RSV immunisation doses per scenario.
# Doses are derived from projected births × scenario uptake and are
# zero outside the vaccination windows. Use this to sanity-check that
# the window logic and birth projections look sensible before submitting.
plot_dose_schedule <- function(doses_df) {
  doses_df %>%
    distinct(scenario_id, target_end_date, value) %>%
    ggplot(aes(x = target_end_date, y = value, colour = scenario_id)) +
    geom_line(linewidth = 0.8) +
    labs(
      title    = "Weekly administered RSV doses by scenario",
      subtitle = "Projected births × scenario uptake, zeroed outside vaccination windows",
      x = "Week ending", y = "Doses", colour = "Scenario"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}
