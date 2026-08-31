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
# scenario:  a scenario_id from the config's `scenarios` block
# age_order: display order for the age bands, normally cfg$age_group_order.
#            Defaults to alphabetical, which orders age bands wrongly -
#            pass the config value to get a sensible axis.
plot_age_breakdown <- function(submission_df, scenario = "no_vacc",
                               age_order = NULL) {

  # pop_group format is "<age>_immTotal" for age-specific total rows
  bands <- submission_df %>%
    filter(target      == "rsv_hospitalisations",
           scenario_id == scenario,
           grepl("_immTotal$", pop_group),
           !grepl("^total_", pop_group)) %>%
    mutate(age_group = sub("_immTotal$", "", pop_group))

  if (is.null(age_order)) age_order <- sort(unique(bands$age_group))

  bands %>%
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


# ---- Plot 4: Adult campaign coverage and protection ---------------------
# Shows the three quantities produced by the adult convolution, for one
# eligible age band:
#
#   coverage       – cumulative fraction of the eligible population
#                    vaccinated. Rises during each campaign window, then
#                    flat (one-off dosing means it never decays).
#   residual VE    – coverage-weighted mean protection among the
#                    vaccinated. Falls as the cohort ages past its dose.
#   effective      – coverage x residual VE, the term that actually
#                    multiplies admissions.
#
# Use this to check that campaign windows, accrual and waning line up:
# coverage should be a staircase over the campaign, residual VE should
# start near the month-0 VE and decay, and effective protection should
# peak shortly after the campaign ends.
#
# Band shows the 95 % interval across the waning-curve ensemble.
plot_adult_protection <- function(protection_df, age_group_to_plot = NULL) {

  d <- protection_df
  if (is.null(age_group_to_plot)) age_group_to_plot <- d$age_group[1]
  d <- d %>% filter(age_group == age_group_to_plot)

  d %>%
    mutate(effective = coverage * residual_ve) %>%
    tidyr::pivot_longer(c("coverage", "residual_ve", "effective"),
                        names_to = "quantity", values_to = "v") %>%
    mutate(quantity = factor(quantity,
                             levels = c("coverage", "residual_ve", "effective"),
                             labels = c("Coverage", "Residual VE",
                                        "Effective (coverage x VE)"))) %>%
    group_by(quantity, target_end_date) %>%
    summarise(median = median(v),
              lo     = quantile(v, 0.025),
              hi     = quantile(v, 0.975),
              .groups = "drop") %>%
    ggplot(aes(x = target_end_date, colour = quantity, fill = quantity)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.8) +
    labs(
      title    = paste0("Adult programme — age group: '", age_group_to_plot, "'"),
      subtitle = "Median + 95 % interval across the waning-curve ensemble",
      x = "Week ending", y = "Fraction", colour = NULL, fill = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}


# ---- Plot 5: Weekly administered dose schedule --------------------------
# Shows the weekly number of RSV immunisation doses per scenario, summed
# across both programmes:
#   infant – births x that scenario's infant uptake, masked to the
#            vaccination windows
#   adult  – eligible population x that week's new coverage increment
#
# For the 2026/27 round infant uptake is zero, so the whole series is the
# adult campaign: a flat block across the campaign weeks and zero either
# side. Use this to check that the campaign window and the population
# denominator line up.
# by_age = FALSE plots the all-ages "undefined" row; TRUE breaks it down
# by the age band the doses were given to.
plot_dose_schedule <- function(doses_df, by_age = FALSE) {

  if (by_age) {
    return(
      doses_df %>%
        filter(pop_group != "undefined") %>%
        distinct(scenario_id, target_end_date, pop_group, value) %>%
        ggplot(aes(x = target_end_date, y = value, colour = pop_group)) +
        geom_line(linewidth = 0.8) +
        facet_wrap(~ scenario_id, ncol = 1, scales = "free_y") +
        labs(title    = "Weekly administered RSV doses by age group",
             subtitle = "Eligible population of each band x that week's coverage increment",
             x = "Week ending", y = "Doses", colour = "Age group") +
        theme_bw() + theme(legend.position = "bottom")
    )
  }

  doses_df %>%
    filter(pop_group == "undefined") %>%
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
