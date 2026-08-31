# ============================================================ #
# Scenario impact plots
#
# One plot per impact table, faceted by scenario with countries on the x
# axis: median point and 90% interval across paired samples, matching
# the ODE model's presentation so the two can be read side by side.
#
# Each takes the corresponding element of compute_scenario_impact().
# ============================================================ #


# Shared skeleton: median + interval, faceted by scenario, coloured by
# season. Season is a single level in this model, so the legend is
# dropped unless there is genuinely more than one.
impact_plot_base <- function(df, title, subtitle, y_lab,
                             zero_line = TRUE, log_y = FALSE) {

  multi_season <- length(unique(df$season)) > 1
  df <- df %>% mutate(season = ifelse(season == "", "all", season))

  p <- ggplot(df, aes(x = iso, y = median,
                      colour = if (multi_season) season else NULL)) +
    { if (zero_line)
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") } +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   position = position_dodge(width = 0.5), linewidth = 0.6) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    facet_wrap(~ scen, ncol = 1, scales = "free_y") +
    labs(title = title, subtitle = subtitle, x = NULL, y = y_lab,
         colour = "Season") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = if (multi_season) "bottom" else "none",
          strip.background = element_rect(fill = "grey90"))

  if (log_y) {
    p <- p + scale_y_log10(labels = scales::comma)
  } else {
    p <- p + scale_y_continuous(labels = scales::comma)
  }
  p
}


# ---- Plot: absolute hospitalisations averted -----------------------------
# Log y axis, because national totals span three orders of magnitude
# between Malta and Germany. Countries with a zero-width interval have
# only one distinct paired difference.
plot_impact_abs_averted <- function(tbl) {
  impact_plot_base(
    tbl %>% filter(median > 0),
    "Hospitalisations averted vs baseline (scenario's eligible ages)",
    "Median and 90% interval across paired samples  -  log10 y axis",
    "Hospitalisations averted (count)",
    zero_line = FALSE, log_y = TRUE)
}


# ---- Plot: relative change vs baseline ----------------------------------
# Negative means averted.
#
# The shaded band spans the waning range over the first year: its lower
# edge is coverage x VE(month 0), a fresh dose; its upper edge is
# coverage x VE(month 12), a dose a year old; the dashed line is the
# midpoint. A real campaign lands between the two, because by the time
# any given admission occurs its doses span a range of ages.
#
# The band's width is WANING; the per-country error bars are Monte-Carlo
# uncertainty. They mean different things and should not be read as one.
plot_impact_pct_averted <- function(tbl, expected = NULL,
                                    label = "scenario's eligible ages") {

  p <- impact_plot_base(
    tbl,
    paste0("Relative change in seasonal RSV hospitalisations vs baseline (",
           label, ")"),
    paste0("Median and 90% interval across paired samples",
           if (!is.null(expected))
             "  -  band: coverage x mean VE, month 0 (lower) to month 12 (upper)"
           else ""),
    "Relative change (%)  -  negative = averted")

  # Reference band drawn UNDER the points, so it never hides them.
  if (!is.null(expected) && nrow(expected) > 0) {
    exp_df <- expected %>% filter(scen %in% unique(tbl$scen))
    p$layers <- c(
      geom_rect(data = exp_df,
                aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
                inherit.aes = FALSE, alpha = 0.18, fill = "steelblue"),
      geom_hline(data = exp_df, aes(yintercept = expected),
                 linetype = "dashed", colour = "steelblue4"),
      p$layers)
  }
  p
}


# ---- Plot: doses per 100k ------------------------------------------------
# Per TOTAL population this varies by country, because the share of the
# population that is 65+ varies. Per ELIGIBLE population it should be
# almost flat at 100000 x coverage - a country far off that line means
# its eligible denominator disagrees with the coverage actually applied.
plot_impact_doses <- function(tbl, per = c("total", "eligible")) {
  per <- match.arg(per)
  impact_plot_base(
    tbl,
    paste0("Vaccine doses administered, per 100k ", per, " population"),
    "Median and 90% interval across paired samples",
    paste0("Doses per 100k ", per, " population"))
}


# ---- Plot: averted per 100k ---------------------------------------------
# Per total population this is the national public-health return; per
# eligible population it is the return among those actually vaccinated,
# which strips out how large each country's 65+ share is.
plot_impact_averted_per100k <- function(tbl, per = c("total", "eligible", "union")) {
  per <- match.arg(per)
  lab <- switch(per,
                total    = "per 100k total population",
                eligible = "per 100k eligible population (scenario's eligible ages)",
                union    = "per 100k eligible population (union of eligible ages)")
  impact_plot_base(
    tbl,
    paste0("Hospitalisations averted ", lab),
    "Median and 90% interval across paired samples",
    paste0("Averted per 100k ",
           if (per == "total") "total" else "eligible", " population"))
}


# Render every impact plot and optionally save them.
#
# Returns a named list of ggplot objects.
build_impact_plots <- function(impact, expected = NULL, dir = NULL,
                               width = 12, height = 10) {

  p <- list(
    abs_averted        = plot_impact_abs_averted(impact$scenario_impact_3_abs_averted),
    pct_scenario_ages  = plot_impact_pct_averted(impact$scenario_impact_1_pct_averted_scenario_ages,
                                                 expected, "scenario's eligible ages"),
    pct_union_ages     = plot_impact_pct_averted(impact$scenario_impact_2_pct_averted_union_ages,
                                                 expected, "union of eligible ages"),
    doses_per100k_tot  = plot_impact_doses(impact$scenario_impact_4a_doses_per100k_total, "total"),
    doses_per100k_elig = plot_impact_doses(impact$scenario_impact_4b_doses_per100k_eligible, "eligible"),
    averted_per100k_tot   = plot_impact_averted_per100k(impact$scenario_impact_5_averted_per100k_total, "total"),
    averted_per100k_elig  = plot_impact_averted_per100k(impact$scenario_impact_6_averted_per100k_eligible, "eligible"),
    averted_per100k_union = plot_impact_averted_per100k(impact$scenario_impact_7_averted_per100k_union, "union")
  )

  if (!is.null(dir)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    for (nm in names(p)) {
      ggsave(file.path(dir, paste0("scenario_impact_", nm, ".png")),
             p[[nm]], width = width, height = height, dpi = 150)
    }
    message("Wrote ", length(p), " impact plots to ", dir)
  }
  p
}
