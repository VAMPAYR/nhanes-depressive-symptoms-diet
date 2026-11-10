# Reporting descriptives generated after the completed paired outcome analysis.
# No new question, outcome fit, model selection, or energy reanalysis.
pm_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
pm_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
pm_equal <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-9))

pm_blank <- function(domain, domain_label, characteristic, source_variable, category, statistic) data.frame(
  domain = domain, domain_label = domain_label, characteristic = characteristic,
  source_variable = source_variable, category = category, statistic = statistic,
  unweighted_n = NA_integer_, domain_n = NA_integer_, denominator_n = NA_integer_, missing_n = NA_integer_,
  weight_denominator = NA_real_, full_domain_weight = NA_real_, weighted_mean = NA_real_, weighted_SD = NA_real_,
  mean_SE = NA_real_, mean_CI_lower = NA_real_, mean_CI_upper = NA_real_,
  weighted_proportion = NA_real_, weighted_percent = NA_real_, proportion_SE = NA_real_,
  proportion_CI_lower = NA_real_, proportion_CI_upper = NA_real_,
  percent_CI_lower = NA_real_, percent_CI_upper = NA_real_, minimum = NA_real_, maximum = NA_real_,
  CI_status = "", domain_PSUs = NA_integer_, domain_strata = NA_integer_, domain_df = NA_real_,
  category_PSUs = NA_integer_, category_strata = NA_integer_, parent_n = NA_integer_,
  parent_PSUs = NA_integer_, parent_strata = NA_integer_, parent_df = NA_real_, stringsAsFactors = FALSE)

pm_add_support <- function(row, support, full_support) {
  row$domain_n <- full_support$n; row$denominator_n <- support$n
  row$missing_n <- full_support$n - support$n
  row$weight_denominator <- support$weight_denominator; row$full_domain_weight <- full_support$weight_denominator
  row$domain_PSUs <- support$psus; row$domain_strata <- support$strata; row$domain_df <- support$domain_df
  row$parent_n <- support$parent_n; row$parent_PSUs <- support$parent_psus
  row$parent_strata <- support$parent_strata; row$parent_df <- support$parent_df
  row
}

pm_continuous <- function(design, mask, variable, characteristic, domain, label, paired) {
  full <- paired$pp_support(design, mask)
  value <- design$variables[[variable]]
  selected <- mask & is.finite(value)
  working <- design; working$variables$.pm_value <- value
  summary <- paired$pp_joint_mean(working, selected, ".pm_value")$summary
  row <- pm_blank(domain, label, characteristic, variable, "", "Mean; weighted SD")
  row <- pm_add_support(row, summary, full)
  row$unweighted_n <- summary$n; row$weighted_mean <- summary$estimate
  row$mean_SE <- summary$se; row$mean_CI_lower <- summary$ci_lower; row$mean_CI_upper <- summary$ci_upper
  row$CI_status <- summary$ci_status
  if (summary$n > 1L) {
    variance <- survey::svyvar(~.pm_value, working[selected, ], na.rm = FALSE)
    row$weighted_SD <- sqrt(as.numeric(stats::coef(variance)))
  }
  if (any(selected)) { row$minimum <- min(value[selected]); row$maximum <- max(value[selected]) }
  pm_assert(pm_equal(row$weighted_mean, weighted.mean(value[selected], design$variables$.weight[selected])),
    paste("Survey and direct weighted means disagree:", domain, variable))
  row
}

pm_categories <- function(design, mask, variable, labels, characteristic, domain, label, paired) {
  value <- as.character(design$variables[[variable]])
  available <- !is.na(value) & value %in% names(labels)
  selected <- mask & available
  full <- paired$pp_support(design, mask)
  table <- do.call(rbind, lapply(names(labels), function(level) {
    summary <- paired$pp_proportion(design, selected, value %in% level)
    row <- pm_blank(domain, label, characteristic, variable, unname(labels[level]), "n; weighted percent")
    row <- pm_add_support(row, summary, full)
    row$unweighted_n <- summary$category_n; row$weighted_proportion <- summary$estimate
    row$weighted_percent <- 100 * summary$estimate; row$proportion_SE <- summary$se
    row$proportion_CI_lower <- summary$ci_lower; row$proportion_CI_upper <- summary$ci_upper
    row$percent_CI_lower <- 100 * summary$ci_lower; row$percent_CI_upper <- 100 * summary$ci_upper
    row$CI_status <- summary$ci_status; row$category_PSUs <- summary$category_psus; row$category_strata <- summary$category_strata
    row
  }))
  pm_assert(sum(table$unweighted_n) == sum(selected) && pm_equal(sum(table$weighted_proportion), 1),
    paste("Categorical counts or weighted proportions fail to reconcile:", domain, variable))
  table
}

pm_table <- function(design, mask, domain, label, paired, full_table = TRUE) {
  support <- paired$pp_support(design, mask)
  participants <- pm_blank(domain, label, "Participants", "PERSON_KEY", "", "Count")
  participants <- pm_add_support(participants, support, support)
  participants$unweighted_n <- support$n; participants$CI_status <- "counts_only"
  rows <- list(participants,
    pm_continuous(design, mask, "RIDAGEYR", "Age, years (80 represents age 80 or older)", domain, label, paired),
    pm_categories(design, mask, "pm_age_band", setNames(c("20-39", "40-59", "60-79", "80+"), c("20-39", "40-59", "60-79", "80+")),
      "Age group, years", domain, label, paired),
    pm_categories(design, mask, "RIAGENDR", c("1" = "Male", "2" = "Female"), "Sex", domain, label, paired))
  if (full_table) rows <- c(rows, list(
    pm_categories(design, mask, "RELEASE_ID", setNames(names(paired$pp_cycles()), names(paired$pp_cycles())), "NHANES release", domain, label, paired),
    pm_continuous(design, mask, "PHQ9_TOTAL", "PHQ-9 symptom score", domain, label, paired),
    pm_categories(design, mask, "pd_phq_band", setNames(c("0-4", "5-9", ">=10"), c("0-4", "5-9", ">=10")),
      "PHQ-9 score group (descriptive)", domain, label, paired),
    pm_categories(design, mask, "pd_order_sign", c("-1" = "First recall day (MEC)", "1" = "Second recall day (telephone)"),
      "Beverage-positive recall day", domain, label, paired),
    pm_continuous(design, mask, "pd_gap_days", "Observed intake-day interval, days", domain, label, paired),
    pm_categories(design, mask, "pm_gap_band", c("below3" = "Below 3 days", "days3to10" = "3-10 days",
      "above10" = "Above 10 days", "unavailable" = "Unreported interval"), "Intake-day interval", domain, label, paired)))
  do.call(rbind, rows)
}

