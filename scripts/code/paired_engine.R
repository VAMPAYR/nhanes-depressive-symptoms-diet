# Source-only paired-recall estimator. No data reads, writes, or fits occur on load.

pe_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

pd_stage_specs <- function() data.frame(
  stage_id = c("primary", "all_intervals", "quadratic"),
  domain_id = c("primary_3_10", "all_intervals", "primary_3_10"),
  quadratic_PHQ = c(FALSE, FALSE, TRUE),
  primary_test = c(TRUE, FALSE, FALSE), stringsAsFactors = FALSE)

pe_config <- function(config) {
  pe_assert(is.list(config) && identical(config$residual_df_policy, "survey_default"),
    "The fixed inference policy requires the actual survey::svyglm residual degrees of freedom.")
  invisible(TRUE)
}

pe_core <- function(core) {
  needed <- c("parent_design", "standardized_contrasts", "design_matrix", "parent_score_covariance")
  pe_assert((is.environment(core) || is.list(core)) &&
    all(vapply(needed, function(name) is.function(core[[name]]), logical(1))),
    "The source-only original survey estimator helpers are required.")
  invisible(TRUE)
}

pe_prepare <- function(master) {
  required <- c("PERSON_KEY", "PHQ9_TOTAL", "RIDAGEYR", "RIAGENDR", "RELEASE_ID", "WTDR2D",
    "SDMVSTRA", "SDMVPSU", "pd_base_eligible", "pd_discordant", "pd_gap_days",
    "pd_delta_other_kcal", "pd_order_sign", paste0("pd_weekday", 2:7))
  pe_assert(is.data.frame(master) && nrow(master) > 0 && all(required %in% names(master)),
    "The full paired-day master lacks required fields.")
  pe_assert(!anyNA(master$PERSON_KEY) && !anyDuplicated(master$PERSON_KEY),
    "Each master row must identify one unique person and one paired-recall record.")
  pe_assert(is.logical(master$pd_base_eligible) && is.logical(master$pd_discordant) &&
    !anyNA(master$pd_base_eligible) && !anyNA(master$pd_discordant),
    "The builder must provide complete logical base-eligibility and discordance flags.")
  numeric_fields <- c("PHQ9_TOTAL", "RIDAGEYR", "RIAGENDR", "WTDR2D", "SDMVSTRA", "SDMVPSU",
    "pd_gap_days", "pd_delta_other_kcal", "pd_order_sign", paste0("pd_weekday", 2:7))
  pe_assert(all(vapply(master[numeric_fields], is.numeric, logical(1))),
    "Canonical model and design fields must use numeric source types.")
  data <- master
  data$.row <- seq_len(nrow(data))
  data$pe_phq5 <- data$PHQ9_TOTAL / 5
  data$pe_phq5_squared <- data$pe_phq5^2
  data$pe_age10 <- (data$RIDAGEYR - 50) / 10
  data$pe_age10_squared <- data$pe_age10^2
  data$pe_sex <- factor(data$RIAGENDR, levels = c(1, 2))
  data$pe_cycle <- factor(data$RELEASE_ID,
    levels = c("2011-2012", "2013-2014", "2015-2016", "2017-2018"))
  parent <- is.finite(data$WTDR2D) & data$WTDR2D > 0 &
    is.finite(data$SDMVSTRA) & is.finite(data$SDMVPSU)
  fields <- c("pe_phq5", "pe_phq5_squared", "pe_age10", "pe_age10_squared", "pe_sex", "pe_cycle",
    "pd_order_sign", paste0("pd_weekday", 2:7), "pd_delta_other_kcal")
  valid <- stats::complete.cases(data[, fields, drop = FALSE])
  for (name in fields) if (is.numeric(data[[name]])) valid <- valid & is.finite(data[[name]])
  weekdays <- as.matrix(data[paste0("pd_weekday", 2:7)])
  valid_weekdays <- rowSums(is.finite(weekdays) & weekdays %in% c(-1, 0, 1)) == 6L &
    rowSums(weekdays == 1, na.rm = TRUE) <= 1L & rowSums(weekdays == -1, na.rm = TRUE) <= 1L
  valid <- valid & data$pd_order_sign %in% c(-1, 1) & valid_weekdays &
    is.finite(data$RIDAGEYR) & data$RIDAGEYR >= 20 &
    is.finite(data$PHQ9_TOTAL) & data$PHQ9_TOTAL >= 0 & data$PHQ9_TOTAL <= 27 &
    data$PHQ9_TOTAL == floor(data$PHQ9_TOTAL)
  valid[is.na(valid)] <- FALSE
  all_intervals <- parent & data$pd_base_eligible & data$pd_discordant & valid
  gap <- is.finite(data$pd_gap_days) & data$pd_gap_days >= 3 & data$pd_gap_days <= 10
  primary <- all_intervals & gap
  masks <- list(primary_3_10 = primary, all_intervals = all_intervals)
  counts <- data.frame(eligibility = c("master_rows", "positive_WTDR2D_parent", "builder_base_eligible",
    "base_eligible_discordant", "all_intervals", "primary_3_10", "all_intervals_missing_gap",
    "all_intervals_observed_gap_outside_3_10", "discordant_excluded_invalid_model_or_design_fields"),
    n = c(nrow(data), sum(parent), sum(data$pd_base_eligible), sum(data$pd_base_eligible & data$pd_discordant),
      sum(all_intervals), sum(primary), sum(all_intervals & !is.finite(data$pd_gap_days)),
      sum(all_intervals & is.finite(data$pd_gap_days) & !gap),
      sum(data$pd_base_eligible & data$pd_discordant & !(parent & valid))), stringsAsFactors = FALSE)
  list(data = data, masks = masks, counts = counts)
}

pd_domains <- function(master) {
  prepared <- pe_prepare(master)
  prepared[c("masks", "counts")]
}

pe_covariance <- function(covariance, label) {
  pe_assert(is.matrix(covariance) && nrow(covariance) == ncol(covariance) && all(is.finite(covariance)),
    paste(label, "must be a finite square matrix."))
  tolerance <- 1e-8 * max(1, max(abs(covariance)))
  # Profile means can be algebraically dependent; positive semidefiniteness is sufficient.
  pe_assert(max(abs(covariance - t(covariance))) <= tolerance &&
    min(eigen((covariance + t(covariance)) / 2, symmetric = TRUE, only.values = TRUE)$values) >= -tolerance,
    paste(label, "is asymmetric or has a negative variance direction."))
  invisible(TRUE)
}

pe_fit_object <- function(master, core, stage_id, config = list(residual_df_policy = "survey_default")) {
  pe_core(core)
  pe_config(config)
  specs <- pd_stage_specs()
  pe_assert(is.character(stage_id) && length(stage_id) == 1L && stage_id %in% specs$stage_id,
    "Only primary, all_intervals, and quadratic are fixed outcome-model stages.")
  spec <- specs[match(stage_id, specs$stage_id), , drop = FALSE]
  prepared <- pe_prepare(master)
  data <- prepared$data
  data$.pe_domain <- prepared$masks[[spec$domain_id]]
  pe_assert(any(data$.pe_domain), "The fixed paired-recall model domain is empty.")
  predictors <- c("pe_phq5", if (spec$quadratic_PHQ) "pe_phq5_squared", "pe_sex", "pe_age10",
    "pe_age10_squared", "pe_cycle", paste0("pd_weekday", 2:7), "pd_order_sign")
  formula <- stats::reformulate(predictors, response = "pd_delta_other_kcal")
  parent <- core$parent_design(data, weight = "WTDR2D")
  sub <- subset(parent, .pe_domain)
  sample <- sub$variables
  frame <- stats::model.frame(formula, sample, na.action = stats::na.fail)
  X <- stats::model.matrix(formula, frame)
  pe_assert(all(is.finite(X)), "The fixed paired-recall model matrix contains nonfinite values.")
  rank <- qr(X, tol = 1e-10)$rank
  pe_assert(rank == ncol(X), "The fixed paired-recall model lost rank. No collapsing or simpler-model fallback is permitted.")
  warnings <- character()
  fit <- withCallingHandlers(survey::svyglm(formula, design = sub, family = stats::gaussian(link = "identity"),
    na.action = stats::na.fail, control = stats::glm.control(epsilon = 1e-10, maxit = 100)),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    })
  pe_assert(isTRUE(fit$converged) && fit$rank == rank &&
    identical(names(stats::coef(fit)), colnames(X)) && all(is.finite(stats::coef(fit))),
    "The fixed paired-recall fit failed to converge or changed its coefficient space.")
  pe_covariance(stats::vcov(fit), "Coefficient covariance")
  residual_df <- fit$df.residual
  pe_assert(is.numeric(residual_df) && length(residual_df) == 1L && is.finite(residual_df) && residual_df > 0,
    "The fitted survey model has nonpositive or invalid stored residual degrees of freedom.")
  diagnostic <- list(stage_id = stage_id, domain_id = spec$domain_id,
    parent_n = nrow(parent$variables), domain_n = nrow(sample),
    parent_strata = length(unique(parent$variables$.stratum)),
    parent_psus = nrow(unique(parent$variables[c(".stratum", ".psu")])),
    parent_df = survey::degf(parent), subdomain_df = survey::degf(sub),
    rank = rank, p = ncol(X), residual_df = residual_df,
    residual_df_policy = "Actual survey::svyglm fit$df.residual; no replacement",
    model_columns = colnames(X), formula = paste(deparse(formula), collapse = " "),
    formula_term_labels = attr(stats::terms(formula), "term.labels"),
    converged = isTRUE(fit$converged), fit_iterations = fit$iter, warnings = unique(warnings),
    weight = "WTDR2D / 4", outcome_scale = "kcal: beverage-positive day minus exact beverage-zero day",
    primary_predictor_scale = "PHQ points / 5")
  list(object = list(fit = fit, data = sample, parent_design = parent, weight = "WTDR2D", diagnostic = diagnostic),
    stage_id = stage_id, spec = spec, domain_counts = prepared$counts)
}

pe_residual_diagnostics <- function(object, core) {
  fitted <- as.numeric(core$design_matrix(object$fit, object$data) %*% stats::coef(object$fit))
  observed <- object$data$pd_delta_other_kcal
  weight <- object$data$.weight
  residual <- observed - fitted
  pe_assert(all(is.finite(fitted)) && all(is.finite(observed)) && all(is.finite(weight) & weight > 0),
    "Residual diagnostics contain nonfinite fitted values, outcomes, or weights.")
  summarize <- function(rows) {
    total <- sum(weight[rows])
    mean_at <- function(value) if (total > 0) sum(weight[rows] * value[rows]) / total else NA_real_
    data.frame(n = sum(rows), weighted_fraction = total / sum(weight),
      observed_mean = mean_at(observed), fitted_mean = mean_at(fitted), residual_bias = mean_at(residual),
      weighted_RMSE = sqrt(mean_at(residual^2)), weighted_MAE = mean_at(abs(residual)), row.names = NULL)
  }
  groups <- list(sex = object$data$pe_sex, cycle = object$data$pe_cycle,
    recall_order_sign = factor(object$data$pd_order_sign, levels = c(-1, 1)),
    PHQ_group = cut(object$data$PHQ9_TOTAL, breaks = c(-Inf, 4, 9, Inf),
      labels = c("0-4", "5-9", ">=10")))
  tables <- list()
  for (grouping in names(groups)) for (level in levels(groups[[grouping]])) {
    tables[[length(tables) + 1L]] <- cbind(data.frame(grouping = grouping, group = level),
      summarize(as.character(groups[[grouping]]) == level))
  }
  # For this Gaussian identity fit, the fitted weighted-QR hat diagonal is
  # diag(sqrt(W) X (X' W X)^-1 X' sqrt(W)); weight rescaling leaves it unchanged.
  hat <- as.numeric(stats::hatvalues(object$fit))
  rank <- object$diagnostic$rank
  pe_assert(length(hat) == length(weight) && all(is.finite(hat)) &&
    all(hat >= -1e-8 & hat <= 1 + 1e-8) &&
    abs(sum(hat) - rank) <= 1e-7 * max(1, rank),
    "Weighted WLS leverage is nonfinite, outside [0,1], or does not sum to the model rank.")
  threshold <- 2 * rank / length(hat)
  above <- hat > threshold
  quantiles <- as.numeric(stats::quantile(hat, probs = c(.25, .5, .75, .9, .95, .99), names = FALSE))
  leverage <- data.frame(n = length(hat), rank = rank, hat_sum = sum(hat),
    hat_min = min(hat), hat_max = max(hat), hat_q25 = quantiles[1], hat_median = quantiles[2],
    hat_q75 = quantiles[3], hat_q90 = quantiles[4], hat_q95 = quantiles[5], hat_q99 = quantiles[6],
    threshold_2p_over_n = threshold, n_above_threshold = sum(above),
    fraction_above_threshold = mean(above), weight_sum = sum(weight),
    weight_sum_above_threshold = sum(weight[above]),
    weighted_fraction_above_threshold = sum(weight[above]) / sum(weight),
    method = "Weighted WLS hat diagonal from the existing Gaussian fit; n is fitted rows and p is model rank",
    threshold_policy = "h > 2p/n is a diagnostic flag only; no deletion, refit, or model change",
    weight_mass_units = "Sum of WTDR2D / 4", row.names = NULL)
  list(description = "Descriptive in-sample residual summaries, with no validated-prediction claim.",
    overall = summarize(rep(TRUE, length(fitted))), grouped = do.call(rbind, tables),
    leverage = leverage,
    fitted_min = min(fitted), fitted_max = max(fitted), residual_min = min(residual), residual_max = max(residual))
}

pe_variance_contributions <- function(influence, object, covariance, estimand_ids) {
  parent <- object$parent_design$variables
  domain <- object$data
  pe_assert(is.matrix(influence) && nrow(influence) == nrow(domain) &&
    ncol(influence) == length(estimand_ids) && all(is.finite(influence)),
    "Every reported estimand requires finite domain influence.")
  pe_assert(!anyDuplicated(parent$.row) && !anyDuplicated(domain$.row), "Duplicate row identities in the influence carrier.")
  positions <- match(domain$.row, parent$.row)
  pe_assert(!anyNA(positions), "A domain influence row is absent from the full WTDR2D parent.")
  carrier <- matrix(0, nrow(parent), ncol(influence))
  carrier[positions, ] <- influence
  psu_tables <- list()
  stratum_tables <- list()
  for (stratum in unique(parent$.stratum)) {
    rows <- which(parent$.stratum == stratum)
    totals <- rowsum(carrier[rows, , drop = FALSE], parent$.psu[rows], reorder = TRUE)
    m <- nrow(totals)
    pe_assert(m > 1L, "The full WTDR2D parent contains a singleton-PSU stratum.")
    centered <- sweep(totals, 2, colMeans(totals))
    contribution <- m / (m - 1) * centered^2
    for (j in seq_along(estimand_ids)) {
      psu_tables[[length(psu_tables) + 1L]] <- data.frame(estimand_id = estimand_ids[j],
        stratum = as.character(stratum), psu = rownames(totals), variance_contribution = contribution[, j], row.names = NULL)
      stratum_tables[[length(stratum_tables) + 1L]] <- data.frame(estimand_id = estimand_ids[j],
        stratum = as.character(stratum), psus = m, variance_contribution = sum(contribution[, j]), row.names = NULL)
    }
  }
  psus <- do.call(rbind, psu_tables)
  strata <- do.call(rbind, stratum_tables)
  psus$variance_share <- NA_real_
  strata$variance_share <- NA_real_
  summaries <- list()
  top_tables <- list()
  for (j in seq_along(estimand_ids)) {
    pr <- which(psus$estimand_id == estimand_ids[j])
    sr <- which(strata$estimand_id == estimand_ids[j])
    total <- sum(psus$variance_contribution[pr])
    saved <- covariance[j, j]
    tolerance <- 1e-12 + 1e-8 * max(abs(total), abs(saved))
    pe_assert(is.finite(total) && total >= 0 && is.finite(saved) &&
      abs(total - saved) <= tolerance && abs(sum(strata$variance_contribution[sr]) - total) <= tolerance,
      "Full-parent influence contributions do not sum to the saved estimand variance.")
    positive <- total > 0
    if (positive) {
      psus$variance_share[pr] <- psus$variance_contribution[pr] / total
      strata$variance_share[sr] <- strata$variance_contribution[sr] / total
      pe_assert(abs(sum(psus$variance_share[pr]) - 1) <= 1e-8 &&
        abs(sum(strata$variance_share[sr]) - 1) <= 1e-8,
        "Positive variance contributions do not define shares summing to one.")
    }
    rp <- pr[order(-psus$variance_contribution[pr], psus$stratum[pr], psus$psu[pr])]
    rs <- sr[order(-strata$variance_contribution[sr], strata$stratum[sr])]
    top <- psus[head(rp, 5L), , drop = FALSE]
    top$rank <- if (positive) seq_len(nrow(top)) else NA_integer_
    top_tables[[j]] <- top[, c("estimand_id", "rank", "stratum", "psu", "variance_contribution", "variance_share")]
    summaries[[j]] <- data.frame(estimand_id = estimand_ids[j],
      status = if (positive) "positive_variance" else "zero_variance", saved_variance = saved,
      contribution_sum = total, variance_sum_error = total - saved, parent_psus = length(pr), parent_strata = length(sr),
      maximum_psu_share = if (positive) psus$variance_share[rp[1]] else NA_real_,
      maximum_psu_stratum = if (positive) psus$stratum[rp[1]] else NA_character_,
      maximum_psu_id = if (positive) psus$psu[rp[1]] else NA_character_,
      top_five_psu_share = if (positive) sum(head(psus$variance_share[rp], 5L)) else NA_real_,
      maximum_stratum_share = if (positive) strata$variance_share[rs[1]] else NA_real_,
      maximum_stratum_id = if (positive) strata$stratum[rs[1]] else NA_character_, row.names = NULL)
  }
  list(method = paste("Taylor variance contributions from every original parent PSU, centered within stratum;",
      "m/(m-1) times squared centered total influence; zero influence outside the model domain.",
      "Includes fitted-coefficient and empirical-standardization-distribution uncertainty; no deletion fits."),
    scale = "kcal_squared", summary = do.call(rbind, summaries), strata = strata, top_psus = do.call(rbind, top_tables))
}

pe_interval <- function(estimate, covariance, ids, df) {
  pe_covariance(covariance, "Reported estimand covariance")
  se <- sqrt(pmax(0, diag(covariance)))
  critical <- stats::qt(.975, df)
  data.frame(estimand_id = ids, estimate = as.numeric(estimate), se = se, df = df,
    lower = as.numeric(estimate) - critical * se, upper = as.numeric(estimate) + critical * se,
    scale = "kcal: beverage-positive day minus exact beverage-zero day", row.names = NULL)
}

pe_aggregate_fit <- function(built, core) {
  object <- built$object
  anchors <- c(0, 5, 10)
  profile_ids <- paste0("PHQ", anchors)
  scenarios <- lapply(anchors, function(phq) {
    scenario <- object$data
    scenario$pe_phq5 <- phq / 5
    scenario$pe_phq5_squared <- (phq / 5)^2
    scenario
  })
  names(scenarios) <- profile_ids
  C <- rbind(PHQ5_minus0 = c(-1, 1, 0), PHQ10_minus0 = c(-1, 0, 1))
  colnames(C) <- profile_ids
  combined <- rbind(diag(3L), C)
  ids <- c(paste0("mean:", profile_ids), paste0("contrast:", rownames(C)))
  dimnames(combined) <- list(ids, profile_ids)
  standardized <- core$standardized_contrasts(object, scenarios, combined, return_influence = TRUE)
  pe_assert(all(is.finite(standardized$estimate)) && all(is.finite(standardized$gradients)),
    "A paired-recall standardized estimate or gradient is nonfinite.")
  pe_covariance(standardized$covariance, "Joint profile and contrast covariance")
  profile_covariance <- standardized$covariance[1:3, 1:3, drop = FALSE]
  dimnames(profile_covariance) <- list(profile_ids, profile_ids)
  contrast_covariance <- standardized$covariance[4:5, 4:5, drop = FALSE]
  dimnames(contrast_covariance) <- list(rownames(C), rownames(C))
  beta <- stats::coef(object$fit)
  coefficient_covariance <- stats::vcov(object$fit)
  # All PHQ contrasts are constant across nuisance rows in these fixed models.
  L <- standardized$gradients[4:5, , drop = FALSE]
  tolerance <- 1e-7 * max(1, max(abs(contrast_covariance)))
  pe_assert(max(abs(as.numeric(L %*% beta) - standardized$estimate[4:5])) <= 1e-7 &&
    max(abs(L %*% coefficient_covariance %*% t(L) - contrast_covariance)) <= tolerance,
    "The PHQ contrasts disagree with their joint coefficient contrasts.")
  df <- object$diagnostic$residual_df
  profiles <- pe_interval(standardized$estimate[1:3], profile_covariance, profile_ids, df)
  names(profiles)[1] <- "profile_id"
  profiles$PHQ9_TOTAL <- anchors
  contrasts <- pe_interval(standardized$estimate[4:5], contrast_covariance, rownames(C), df)
  names(contrasts)[1] <- "contrast_id"
  coefficients <- pe_interval(beta, coefficient_covariance, names(beta), df)
  names(coefficients)[1] <- "term"
  coefficients$scale <- "Outcome kcal per model term; PHQ predictor unit is 5 points"
  primary_estimate <- NULL
  primary_wald <- NULL
  if (built$stage_id == "primary") {
    pe_assert(abs(contrasts$estimate[1] - beta["pe_phq5"]) <= 1e-7,
      "The primary five-point profile contrast differs from the PHQ coefficient.")
    variance <- coefficient_covariance["pe_phq5", "pe_phq5"]
    pe_assert(is.finite(variance) && variance > 0, "The primary PHQ coefficient has nonpositive variance.")
    primary_estimate <- coefficients[coefficients$term == "pe_phq5", , drop = FALSE]
    primary_estimate$scale <- "kcal per 5 PHQ points"
    statistic <- as.numeric(beta["pe_phq5"]^2 / variance)
    primary_wald <- list(test = "Primary PHQ trend coefficient equals zero", statistic = statistic,
      numerator_df = 1, denominator_df = df, p_value = stats::pf(statistic, 1, df, lower.tail = FALSE),
      method = "Same-fit survey Wald F; matching t confidence interval")
  }
  list(stage_id = built$stage_id, domain_id = built$spec$domain_id,
    outcome = "pd_delta_other_kcal", outcome_definition = paste("Within-person oriented difference in other food/beverage kcal:",
      "beverage-positive day minus exact beverage-zero day; excludes coded positive-ethanol beverage items and linked beverage-addition energy."),
    domain_counts = built$domain_counts, coefficients = coefficients, coefficient_covariance = coefficient_covariance,
    profiles = profiles, profile_covariance = profile_covariance, contrasts = contrasts, contrast_covariance = contrast_covariance,
    contrast_definitions = C, joint_covariance = standardized$covariance,
    conditional_distribution_covariance = standardized$conditional_distribution_covariance,
    standardized_gradients = standardized$gradients, variance_method = standardized$variance_method,
    standardization_distribution = paste("Own fitted domain:", built$spec$domain_id,
      "; same original empirical nuisance rows and normalized WTDR2D/4 weights at PHQ0,5,10"),
    primary_estimate = primary_estimate, primary_wald = primary_wald, diagnostic = object$diagnostic,
    residual_diagnostics = pe_residual_diagnostics(object, core),
    influence_contributions = pe_variance_contributions(standardized$influence, object, standardized$covariance, ids))
}

pd_fit_stage <- function(master, core, stage_id, config = list(residual_df_policy = "survey_default")) {
  built <- pe_fit_object(master, core, stage_id, config)
  pe_aggregate_fit(built, core)
}
