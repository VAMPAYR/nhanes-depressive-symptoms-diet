# Source-only, fixed secondary analyses. See SECONDARY_SPEC.md.
se_stages <- function() c("PHQ_order", "PHQ_sex", "appetite_deleted")

se_prepare <- function(master) {
  prepared <- pe_prepare(master)
  data <- prepared$data
  data$.se_domain <- prepared$masks$primary_3_10
  assert(sum(data$.se_domain) == 2320L, "The original primary cohort changed.")
  assert(all(data$DPQ050[data$.se_domain] %in% 0:3), "Appetite item missing in the original cohort.")
  data$se_appetite_deleted <- data$PHQ9_TOTAL - data$DPQ050
  data$se_appetite5 <- data$se_appetite_deleted / 5
  assert(all(data$se_appetite_deleted[data$.se_domain] %in% 0:24), "Invalid appetite-deleted score.")
  list(data = data, key_sha256 = digest::digest(data$PERSON_KEY[data$.se_domain], algo = "sha256"),
    counts = prepared$counts)
}

se_fit <- function(data, stage) {
  assert(stage %in% se_stages(), "Unknown secondary stage.")
  nuisance <- c("pe_sex", "pe_age10", "pe_age10_squared", "pe_cycle", paste0("pd_weekday", 2:7), "pd_order_sign")
  predictors <- c(if (stage == "appetite_deleted") "se_appetite5" else "pe_phq5", nuisance,
    if (stage == "PHQ_order") "pe_phq5:pd_order_sign", if (stage == "PHQ_sex") "pe_phq5:pe_sex")
  formula <- stats::reformulate(predictors, response = "pd_delta_other_kcal")
  parent <- parent_design(data, weight = "WTDR2D")
  assert(nrow(parent$variables) == 28845L && survey::degf(parent) == 62L,
    "Original survey parent changed.")
  domain <- subset(parent, .se_domain)
  sample <- domain$variables
  assert(nrow(sample) == 2320L, "Secondary model changed the primary cohort.")
  frame <- stats::model.frame(formula, sample, na.action = stats::na.fail)
  X <- stats::model.matrix(formula, frame)
  assert(all(is.finite(X)) && qr(X, tol = 1e-10)$rank == ncol(X), "Secondary model is nonfinite or rank deficient.")
  warnings <- character()
  fit <- withCallingHandlers(survey::svyglm(formula, design = domain,
    family = stats::gaussian(link = "identity"), na.action = stats::na.fail,
    control = stats::glm.control(epsilon = 1e-10, maxit = 100)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  assert(isTRUE(fit$converged) && fit$rank == ncol(X) &&
    identical(names(stats::coef(fit)), colnames(X)) && all(is.finite(stats::coef(fit))), "Secondary model failed.")
  pe_covariance(stats::vcov(fit), "Secondary coefficient covariance")
  df <- fit$df.residual
  assert(is.finite(df) && df > 0, "Invalid secondary stored residual degrees of freedom.")
  diagnostic <- list(stage_id = stage, domain_n = nrow(sample), parent_n = nrow(parent$variables),
    parent_psus = nrow(unique(parent$variables[c(".stratum", ".psu")])),
    parent_strata = length(unique(parent$variables$.stratum)), parent_df = survey::degf(parent),
    subdomain_df = survey::degf(domain), residual_df = df, rank = fit$rank,
    formula = paste(deparse(formula), collapse = " "), converged = fit$converged,
    warnings = unique(warnings), weight = "WTDR2D / 4", cohort = "Original primary 2320 participants")
  object <- list(fit = fit, data = sample, parent_design = parent, diagnostic = diagnostic)
  L <- matrix(0, if (stage == "appetite_deleted") 2L else 3L, length(stats::coef(fit)),
    dimnames = list(NULL, names(stats::coef(fit))))
  if (stage == "PHQ_order") {
    rownames(L) <- c("first_positive_PHQ10_minus0", "second_positive_PHQ10_minus0", "second_minus_first_PHQ10_contrast")
    L[1:2, "pe_phq5"] <- 2
    L[1:2, "pe_phq5:pd_order_sign"] <- c(-2, 2)
    L[3, ] <- L[2, ] - L[1, ]
    interaction <- "pe_phq5:pd_order_sign"
  } else if (stage == "PHQ_sex") {
    rownames(L) <- c("male_PHQ10_minus0", "female_PHQ10_minus0", "female_minus_male_PHQ10_contrast")
    L[1:2, "pe_phq5"] <- 2
    L[2, "pe_phq5:pe_sex2"] <- 2
    L[3, ] <- L[2, ] - L[1, ]
    interaction <- "pe_phq5:pe_sex2"
  } else {
    rownames(L) <- c("appetite_deleted_per5", "appetite_deleted_10_minus0")
    L[, "se_appetite5"] <- c(1, 2)
  }
  beta <- stats::coef(fit)
  covariance <- L %*% stats::vcov(fit) %*% t(L)
  estimate <- as.numeric(L %*% beta)
  # Each within-group PHQ contrast is constant across nuisance rows. Its
  # empirical-standardization contribution cancels; joint covariance is L V L'.
  scenarios <- list()
  if (stage != "appetite_deleted") {
    for (group in 1:2) for (score in c(0, 10)) {
      z <- sample; z$pe_phq5 <- score / 5
      if (stage == "PHQ_order") z$pd_order_sign <- c(-1, 1)[group]
      else z$pe_sex <- factor(rep(group, nrow(z)), levels = c(1, 2))
      scenarios[[length(scenarios) + 1L]] <- z
    }
    names(scenarios) <- c("group1_PHQ0", "group1_PHQ10", "group2_PHQ0", "group2_PHQ10")
    C <- rbind(c(-1, 1, 0, 0), c(0, 0, -1, 1), c(1, -1, -1, 1))
  } else {
    scenarios <- lapply(c(0, 5, 10), function(score) { z <- sample; z$se_appetite5 <- score / 5; z })
    names(scenarios) <- c("score0", "score5", "score10")
    C <- rbind(c(-1, 1, 0), c(-1, 0, 1))
  }
  rownames(C) <- rownames(L)
  joint <- standardized_contrasts(object, scenarios, C)
  tolerance <- 1e-8 * max(1, max(abs(covariance)))
  assert(max(abs(joint$estimate - estimate)) < 1e-7 &&
    max(abs(joint$covariance - covariance)) < tolerance && max(abs(joint$gradients - L)) < 1e-7,
    "Secondary contrast derivation disagrees with full-parent Taylor influence.")
  table <- pe_interval(estimate, covariance, rownames(L), df)
  table$stage_id <- stage
  table$predictor <- if (stage == "appetite_deleted") "Appetite-deleted score, PHQ9_TOTAL minus DPQ050; no diagnostic thresholds" else "PHQ9_TOTAL"
  coefficients <- pe_interval(beta, stats::vcov(fit), names(beta), df)
  names(coefficients)[1] <- "term"
  coefficients$stage_id <- stage
  test <- NULL
  if (stage != "appetite_deleted") {
    v <- stats::vcov(fit)[interaction, interaction]
    assert(is.finite(v) && v > 0, "Invalid interaction variance.")
    F <- as.numeric(beta[interaction]^2 / v)
    test <- data.frame(stage_id = stage, term = interaction, statistic = F, numerator_df = 1,
      denominator_df = df, p_value = stats::pf(F, 1, df, lower.tail = FALSE),
      method = "One-df survey Wald F; stored svyglm residual df")
  }
  support <- list()
  if (stage != "appetite_deleted") for (group in 1:2) {
    mask <- if (stage == "PHQ_order") parent$variables$.se_domain & parent$variables$pd_order_sign %in% c(-1, 1)[group]
      else parent$variables$.se_domain & parent$variables$RIAGENDR %in% group
    label <- if (stage == "PHQ_order") c("First recall positive", "Second recall positive")[group] else c("Male", "Female")[group]
    for (high in c(FALSE, TRUE)) {
      selected <- mask & if (high) parent$variables$PHQ9_TOTAL >= 10 else TRUE
      selected[is.na(selected)] <- FALSE
      support[[length(support) + 1L]] <- cbind(data.frame(stage_id = stage, group = label,
        PHQ_subset = if (high) ">=10" else "All scores"), pp_support(parent, selected))
    }
  }
  list(table = table, coefficients = coefficients, contrast_covariance = covariance, contrast_matrix = L,
    coefficient_covariance = stats::vcov(fit), test = test, diagnostic = diagnostic,
    support = oa_bind_tables(support), validation = list(full_parent_Taylor_equals_coefficient_contrasts = TRUE,
      maximum_covariance_error = max(abs(joint$covariance - covariance))))
}

se_patterns <- function(data) {
  parent <- parent_design(data, weight = "WTDR2D")
  v <- parent$variables
  patterns <- c("both_zero", "day1_positive_day2_zero", "day1_zero_day2_positive", "both_positive")
  rows <- list()
  for (sex in 1:2) for (band in levels(v$pd_phq_band)) {
    domain <- v$pd_base_eligible & v$RIAGENDR %in% sex & v$pd_phq_band %in% band
    assert(all(v$pd_pair_pattern[domain] %in% patterns), "A descriptive recall pair is unclassified.")
    cell <- list()
    for (pattern in patterns) {
      category <- domain & v$pd_pair_pattern %in% pattern
      z <- pp_support(parent, domain)
      s <- pp_support(parent, category)
      z$category_n <- s$n; z$category_psus <- s$psus; z$category_strata <- s$strata
      z$category_weight <- s$weight_denominator
      z$weighted_proportion <- s$weight_denominator / z$weight_denominator
      z$weighted_percent <- 100 * z$weighted_proportion
      z$missing_pattern_n <- sum(domain & !v$pd_pair_pattern %in% patterns)
      z$sex <- if (sex == 1) "Male" else "Female"
      z$PHQ_band <- band; z$pattern <- pattern
      cell[[length(cell) + 1L]] <- z
    }
    cell <- do.call(rbind, cell)
    assert(sum(cell$category_n) == sum(domain) && abs(sum(cell$weighted_proportion) - 1) < 1e-10,
      "Recall-pattern counts or percentages do not sum to their sex-band denominator.")
    rows[[length(rows) + 1L]] <- cell
  }
  list(table = do.call(rbind, rows), metadata = list(
    domain = "Adults >=20, strict complete PHQ, two reliable valid recalls, positive WTDR2D, valid survey identifiers",
    denominator = "Within each sex and PHQ-band cell", model = "Unadjusted descriptions; no outcome model or risk test",
    intervals = "Counts and weighted percentages only; no confidence intervals or tests"))
}
