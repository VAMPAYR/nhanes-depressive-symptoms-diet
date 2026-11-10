# Reporting orchestration uses unchanged pm_* descriptive functions.
rp_descriptives <- function(data, core, paired) {
  domains <- paired$pd_domains(data)
  z <- data[c("PERSON_KEY", "RIDAGEYR", "RIAGENDR", "RELEASE_ID", "PHQ9_TOTAL", "WTDR2D",
    "SDMVSTRA", "SDMVPSU", "pd_phq_total", "pd_phq_band", "pd_day1_valid", "pd_day2_valid",
    "pd_order_sign", "pd_gap_days")]
  z$.pm_row <- seq_len(nrow(z))
  z$pm_age_band <- as.character(cut(z$RIDAGEYR, breaks = c(19, 39, 59, 79, Inf), labels = c("20-39", "40-59", "60-79", "80+")))
  z$pm_gap_band <- ifelse(!is.finite(z$pd_gap_days), "unavailable", ifelse(z$pd_gap_days < 3, "below3", ifelse(z$pd_gap_days <= 10, "days3to10", "above10")))
  z$.pm_primary <- domains$masks$primary_3_10
  z$.pm_all <- domains$masks$all_intervals
  parent <- core$parent_design(z, weight = "WTDR2D")
  ps <- paired$pp_support(parent, rep(TRUE, nrow(parent$variables)))
  attr(parent, "pp_parent_support") <- list(n = ps$n, psus = ps$psus, strata = ps$strata, df = ps$parent_df)
  v <- parent$variables
  tables <- list(table1_primary = pm_table(parent, v$.pm_primary, "primary_3_10", "Primary: observed interval 3-10 days", paired),
    table1_all_intervals = pm_table(parent, v$.pm_all, "all_intervals", "Sensitivity: all intervals, including unreported", paired))
  adult <- is.finite(v$RIDAGEYR) & v$RIDAGEYR >= 20
  complete <- adult & v$pd_day1_valid & v$pd_day2_valid & is.finite(v$pd_phq_total)
  missing <- adult & v$pd_day1_valid & v$pd_day2_valid & !is.finite(v$pd_phq_total)
  selection <- list(adult_two_day_weight_parent = adult, adult_complete_PHQ = complete, adult_PHQ_unavailable = missing)
  labels <- c(adult_two_day_weight_parent = "Adult positive-two-day-weight parent", adult_complete_PHQ = "Adults with two reliable days and complete strict PHQ",
    adult_PHQ_unavailable = "Adults with two reliable days and unavailable strict PHQ")
  tables$selection_age_sex <- do.call(rbind, lapply(names(selection), function(name) pm_table(parent, selection[[name]], name, labels[[name]], paired, full_table = FALSE)))
  masks <- c(list(full_two_day_weight_parent = rep(TRUE, nrow(v)), primary_3_10 = v$.pm_primary, all_intervals = v$.pm_all), selection)
  tables$domain_support <- do.call(rbind, lapply(names(masks), function(name) cbind(domain = name, paired$pp_support(parent, masks[[name]]))))
  pm_assert(sum(domains$masks$primary_3_10) == 2320L && sum(domains$masks$all_intervals) == 2701L &&
    ps$n == 28845L && ps$psus == 121L && ps$strata == 59L && ps$parent_df == 62 &&
    sum(adult) == 17131L && sum(complete) == 16218L && sum(missing) == 913L && all((complete | missing) == adult),
    "Canonical participant and survey support counts changed.")
  tables
}

rp_publication_figure <- function(d, path) {
  stopifnot(nrow(d) == 9L, setequal(d$stage_id, c("primary", "all_intervals", "quadratic")),
    all(is.finite(as.matrix(d[c("estimate", "lower", "upper")]))))
  grDevices::png(path, width = 1950, height = 1020, res = 300, pointsize = 10,
    type = "cairo", family = "Times New Roman")
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(2.4, 3.5, 2.1, .3), oma = c(1.4, .3, 0, 0),
    mgp = c(2.1, .6, 0), tcl = -.25, cex = 1)
  ids <- c("primary", "all_intervals", "quadratic")
  labels <- c("A  Primary", "B  All intervals", "C  Quadratic")
  for (j in seq_along(ids)) {
    z <- d[d$stage_id == ids[j], ]; z <- z[order(z$PHQ9_TOTAL), ]
    stopifnot(identical(as.numeric(z$PHQ9_TOTAL), c(0, 5, 10)))
    graphics::plot(NA, xlim = c(-.9, 10.9), ylim = c(-150, 250),
      xlab = "", ylab = if (j == 1) "Paired energy difference (kcal)" else "",
      xaxt = "n", yaxt = "n", bty = "l", main = labels[j], cex.main = 1)
    graphics::axis(1, at = c(0, 5, 10))
    graphics::axis(2, at = c(-100, 0, 100, 200), las = 1)
    graphics::abline(h = 0, lty = 2, col = "grey60", lwd = .7)
    graphics::arrows(z$PHQ9_TOTAL, z$lower, z$PHQ9_TOTAL, z$upper,
      angle = 90, code = 3, length = .04, lwd = 1)
    graphics::points(z$PHQ9_TOTAL, z$estimate, pch = 16, cex = .8)
  }
  graphics::mtext("PHQ-9 symptom score", side = 1, outer = TRUE, line = .25)
  invisible(path)
}
