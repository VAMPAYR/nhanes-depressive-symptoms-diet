# Source-only functions. Reading this file does not read data or fit models.
st_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
st_hash <- function(path) digest::digest(file = path, algo = "sha256")
st_key_hash <- function(keys) digest::digest(sort(as.character(keys)), algo = "sha256")
st_near <- function(x, y, tolerance = 1e-7) {
  identical(dim(x), dim(y)) && length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
    max(abs(x - y)) <= tolerance * max(1, max(abs(y)))
}
st_eligible <- function(d) {
  # Rebuild eligibility without using pd_base_eligible, which requires PHQ.
  w <- as.matrix(d[paste0("pd_weekday", 2:7)])
  weekday_ok <- rowSums(is.finite(w) & w %in% c(-1, 0, 1)) == 6L &
    rowSums(w == 1, na.rm = TRUE) <= 1L & rowSums(w == -1, na.rm = TRUE) <= 1L
  z <- is.finite(d$RIDAGEYR) & d$RIDAGEYR >= 20 & is.finite(d$WTDR2D) & d$WTDR2D > 0 &
    is.finite(d$SDMVSTRA) & is.finite(d$SDMVPSU) & d$pd_day1_valid & d$pd_day2_valid &
    d$pd_discordant & is.finite(d$pd_gap_days) & d$pd_gap_days >= 3 & d$pd_gap_days <= 10 &
    d$RIAGENDR %in% c(1, 2) & d$RELEASE_ID %in% names(pp_cycles()) &
    is.finite(d$pd_delta_other_kcal) & d$pd_order_sign %in% c(-1, 1) & weekday_ok
  z[is.na(z)] <- FALSE
  z
}
st_context <- function(d) {
  d$st_income <- ifelse(is.finite(d$INDFMPIR) & d$INDFMPIR >= 0 & d$INDFMPIR <= 5, d$INDFMPIR, NA_real_)
  edu <- ifelse(d$DMDEDUC2 %in% 1:2, "less_than_HS", ifelse(d$DMDEDUC2 == 3, "HS_GED",
    ifelse(d$DMDEDUC2 == 4, "some_college_AA", ifelse(d$DMDEDUC2 == 5, "college_plus", NA))))
  d$st_education <- factor(edu, levels = c("less_than_HS", "HS_GED", "some_college_AA", "college_plus"))
  d$st_race <- factor(ifelse(d$RIDRETH1 %in% 1:5, d$RIDRETH1, NA), levels = 1:5)
  smoking <- ifelse(d$SMQ020 %in% 2, "never", ifelse(d$SMQ020 %in% 1 & d$SMQ040 %in% 3,
    "former", ifelse(d$SMQ020 %in% 1 & d$SMQ040 %in% 1:2, "current", NA)))
  d$st_smoking <- factor(smoking, levels = c("never", "former", "current"))
  d$st_context_complete <- complete.cases(d[c("st_income", "st_education", "st_race", "st_smoking")])
  d
}
st_formula <- function(quadratic = FALSE, context = FALSE) stats::reformulate(c("pe_phq5",
  if (quadratic) "pe_phq5_squared", "pe_sex", "pe_age10", "pe_age10_squared", "pe_cycle",
  paste0("pd_weekday", 2:7), "pd_order_sign",
  if (context) c("st_income", "st_education", "st_race", "st_smoking")), response = "pd_delta_other_kcal")
st_L <- function(columns) {
  L <- setNames(rep(0, length(columns)), columns)
  st_assert("pe_phq5" %in% columns, "PHQ coefficient absent.")
  L["pe_phq5"] <- 2
  if ("pe_phq5_squared" %in% columns) L["pe_phq5_squared"] <- 4
  L
}
st_wls <- function(X, y, w) {
  st_assert(all(is.finite(X)) && all(is.finite(y)) && all(is.finite(w) & w >= 0) && sum(w) > 0,
    "WLS input is invalid.")
  fit <- stats::lm.wfit(X, y, w / mean(w), tol = 1e-10)
  st_assert(fit$rank == ncol(X) && all(is.finite(fit$coefficients)), "Fixed model lost rank.")
  fit
}
st_contrast <- function(beta, V, df, label, n, high, rank) {
  L <- st_L(names(beta)); estimate <- sum(L * beta); variance <- as.numeric(t(L) %*% V %*% L)
  st_assert(is.finite(variance) && variance >= 0 && df > 0, "Invalid contrast inference.")
  se <- sqrt(variance); critical <- stats::qt(.975, df)
  data.frame(model = label, n = n, phq_ge10_n = high, rank = rank, df = df,
    estimate = estimate, SE = se, lower95 = estimate - critical * se,
    upper95 = estimate + critical * se, p_value = 2 * stats::pt(-abs(estimate / se), df))
}
st_delete <- function(X, y, w, stratum, key, original_se, stage) {
  fit <- st_wls(X, y, w); b <- fit$coefficients; L <- st_L(colnames(X))
  wn <- w / mean(w); inverse <- solve(crossprod(X, X * wn)); Q <- X %*% inverse
  h <- rowSums(Q * X) * wn
  st_assert(all(h < 1 - 1e-10), "A case deletion destroys the fixed coefficient space.")
  shift <- -as.numeric(Q %*% L) * wn * fit$residuals / (1 - h)
  full <- sum(L * b)
  cases <- data.frame(stage = stage, person_key = key, stratum = stratum, leverage = h,
    estimate = full + shift, shift = shift, rank_failure = FALSE)
  checked <- unique(c(1L, which.max(h), which.max(abs(shift))))
  for (i in checked) {
    direct <- sum(L * st_wls(X[-i, , drop = FALSE], y[-i], w[-i])$coefficients)
    st_assert(st_near(direct, cases$estimate[i], 1e-8), "Case-deletion identity disagrees with direct refit.")
  }
  strata <- do.call(rbind, lapply(sort(unique(stratum)), function(s) {
    keep <- stratum != s
    f <- tryCatch(st_wls(X[keep, , drop = FALSE], y[keep], w[keep]), error = identity)
    failed <- inherits(f, "error")
    estimate <- if (failed) NA_real_ else sum(L * f$coefficients)
    data.frame(stage = stage, stratum = s, omitted_n = sum(!keep), estimate = estimate,
      shift = estimate - full, rank_failure = failed,
      failure_reason = if (failed) conditionMessage(f) else "")
  }))
  summarize <- function(z, unit) data.frame(stage = stage, deletion_unit = unit, full_estimate = full,
    deletion_min = min(z$estimate, na.rm = TRUE), deletion_max = max(z$estimate, na.rm = TRUE),
    max_absolute_shift = max(abs(z$shift), na.rm = TRUE),
    max_shift_over_original_SE = max(abs(z$shift), na.rm = TRUE) / original_se,
    rank_failures = sum(z$rank_failure), direct_case_refits_checked = if (unit == "person") length(checked) else NA_integer_)
  list(cases = cases, strata = strata, summary = rbind(summarize(cases, "person"), summarize(strata, "stratum")))
}
st_load_data <- function(data_root) {
  paired <- pp_build_data(mb_build(data_root)$data, data_root)
  d <- paired$data; raw <- list()
  d$st_dpq_record_present <- FALSE
  for (cycle in names(pp_cycles())) for (component in c("DEMO", "SMQ", "DPQ")) {
    path <- file.path(data_root, cycle, paste0(component, "_", pp_cycles()[[cycle]], ".xpt"))
    dat <- foreign::read.xport(path)
    st_assert(!anyDuplicated(dat$SEQN), paste("Duplicate source key:", path))
    rows <- which(d$RELEASE_ID == cycle); match_at <- match(d$SEQN[rows], dat$SEQN)
    if (component == "DPQ") d$st_dpq_record_present[rows] <- !is.na(match_at) else {
      fields <- switch(component, DEMO = c("INDFMPIR", "DMDEDUC2", "RIDRETH1"), SMQ = c("SMQ020", "SMQ040"))
      st_assert(all(fields %in% names(dat)), paste("Missing context source field:", path))
      for (field in fields) {
        if (!field %in% names(d)) d[[field]] <- NA_real_
        d[[field]][rows] <- dat[[field]][match_at]
      }
    }
    raw[[length(raw) + 1L]] <- data.frame(path = normalizePath(path, winslash = "/"), sha256 = st_hash(path))
  }
  d <- st_context(pe_prepare(d)$data)
  d$st_eligible <- st_eligible(d)
  items <- as.matrix(d[sprintf("DPQ%03d", seq(10, 90, 10))])
  d$st_valid_items <- rowSums(matrix(items %in% 0:3, nrow = nrow(d)))
  d$st_complete <- d$st_valid_items == 9L
  d$st_all_items_NA <- rowSums(is.na(items)) == 9L
  d$st_any_refused <- rowSums(items == 7, na.rm = TRUE) > 0L
  d$st_any_unknown <- rowSums(items == 9, na.rm = TRUE) > 0L
  st_assert(identical(d$st_eligible & d$st_complete, pe_prepare(d)$masks$primary_3_10),
    "PHQ-independent eligibility fails to recover the original participant keys.")
  list(data = d, raw = do.call(rbind, raw), builder_audit = paired$audit)
}
st_missingness <- function(d) {
  e <- d[d$st_eligible, , drop = FALSE]; e$st_mean_energy <- (e$pd_day1_total_kcal + e$pd_day2_total_kcal) / 2
  e$st_drinking_ethanol <- ifelse(e$pd_order_sign == 1, e$pd_day2_beverage_alcohol_g, e$pd_day1_beverage_alcohol_g)
  continuous <- c("RIDAGEYR", "pd_delta_other_kcal", "st_mean_energy", "st_drinking_ethanol")
  categorical <- c("RIAGENDR", "RELEASE_ID", "pd_order_sign")
  rows <- list(); flow <- list()
  for (complete in c(TRUE, FALSE)) {
    z <- e[e$st_complete == complete, , drop = FALSE]; w <- z$WTDR2D / 4; label <- if (complete) "complete" else "incomplete"
    flow[[label]] <- data.frame(response = label, n = nrow(z), weighted_total = sum(w),
      weighted_fraction = sum(w) / sum(e$WTDR2D / 4), absent_dpq_record_n = sum(!z$st_dpq_record_present),
      present_all_invalid_n = sum(z$st_dpq_record_present & z$st_valid_items == 0),
      all_nine_NA_n = sum(z$st_all_items_NA), partial_n = sum(z$st_valid_items %in% 1:8),
      any_refused_n = sum(z$st_any_refused), any_unknown_n = sum(z$st_any_unknown))
    for (field in continuous) {
      x <- z[[field]]; valid <- is.finite(x); mu <- weighted.mean(x[valid], w[valid])
      rows[[length(rows) + 1L]] <- data.frame(response = label, variable = field, level = "continuous",
        n = sum(valid), weighted_mean_or_proportion = mu,
        weighted_population_SD = sqrt(weighted.mean((x[valid] - mu)^2, w[valid])))
    }
    for (field in categorical) for (level in sort(unique(e[[field]]))) {
      yes <- z[[field]] == level
      rows[[length(rows) + 1L]] <- data.frame(response = label, variable = field, level = as.character(level),
        n = sum(yes), weighted_mean_or_proportion = sum(w[yes]) / sum(w), weighted_population_SD = NA_real_)
    }
  }
  pattern <- aggregate(list(n = rep(1L, nrow(e))), e[c("st_complete", "st_valid_items", "st_dpq_record_present", "st_all_items_NA", "st_any_refused", "st_any_unknown")], sum)
  weighted_pattern <- aggregate(list(weighted_total = e$WTDR2D / 4), e[c("st_complete", "st_valid_items", "st_dpq_record_present", "st_all_items_NA", "st_any_refused", "st_any_unknown")], sum)
  pattern$weighted_total <- weighted_pattern$weighted_total
  pattern$weighted_fraction_of_eligible <- pattern$weighted_total / sum(e$WTDR2D / 4)
  list(flow = do.call(rbind, flow), description = do.call(rbind, rows), item_patterns = pattern)
}
st_fit_context <- function(d, core) {
  domain <- d$st_eligible & d$st_complete & d$st_context_complete
  d$.st_domain <- domain; parent <- core$parent_design(d, "WTDR2D"); sub <- subset(parent, .st_domain)
  fits <- lapply(c(FALSE, TRUE), function(context) {
    formula <- st_formula(context = context); X <- model.matrix(formula, sub$variables)
    st_assert(qr(X, tol = 1e-10)$rank == ncol(X), "Prespecified context matrix lost rank.")
    f <- survey::svyglm(formula, design = sub, family = gaussian(), na.action = na.fail,
      control = glm.control(epsilon = 1e-10, maxit = 100))
    st_assert(f$converged && f$rank == ncol(X), "Context fit failed.")
    list(fit = f, row = st_contrast(coef(f), vcov(f), f$df.residual,
      if (context) "same_sample_extended_context" else "same_sample_original_baseline",
      sum(domain), sum(d$PHQ9_TOTAL[domain] >= 10), f$rank))
  })
  st_assert(identical(rownames(model.frame(fits[[1]]$fit)), rownames(model.frame(fits[[2]]$fit))), "Context models use different cases.")
  list(comparison = do.call(rbind, lapply(fits, `[[`, "row")),
    coefficients = lapply(fits, function(z) coef(z$fit)), covariances = lapply(fits, function(z) vcov(z$fit)),
    participant_key_sha256 = st_key_hash(d$PERSON_KEY[domain]))
}
st_response_fit <- function(R, complete, w) {
  use <- w > 0
  f <- stats::glm.fit(R[use, , drop = FALSE], as.numeric(complete[use]),
    weights = w[use] / mean(w[use]), family = quasibinomial(),
    control = glm.control(epsilon = 1e-10, maxit = 100))
  st_assert(f$converged && f$rank == ncol(R) && all(is.finite(f$coefficients)), "Response model failed convergence or rank.")
  p <- plogis(as.numeric(R %*% f$coefficients))
  st_assert(all(is.finite(p) & p > 0 & p <= 1), "Response probabilities invalid.")
  list(probability = p, coefficients = f$coefficients, rank = f$rank, iterations = f$iter)
}
st_weight_diagnostics <- function(w, label) data.frame(weight = label, n = length(w),
  total = sum(w), minimum = min(w), maximum = max(w), maximum_share = max(w) / sum(w),
  top_one_percent_share = sum(sort(w, decreasing = TRUE)[seq_len(max(1, ceiling(.01 * length(w))))]) / sum(w),
  kish_weight_only_effective_sample_size = sum(w)^2 / sum(w^2))
st_response <- function(d, core, cfg, progress = function(...) invisible(NULL)) {
  parent <- core$parent_design(d, "WTDR2D")
  rep_parent <- survey::as.svrepdesign(parent, type = "JKn", mse = TRUE)
  # Construct replicates on the full parent, then take the eligibility domain.
  rep_eligible <- subset(rep_parent, st_eligible)
  e <- rep_eligible$variables; complete <- e$st_complete
  e$st_outcome1000 <- e$pd_delta_other_kcal / 1000
  response_formula <- ~ pe_age10 + pe_age10_squared + pe_sex + pe_cycle + pd_order_sign + st_outcome1000
  R <- model.matrix(response_formula, e)
  X <- model.matrix(st_formula(), e[complete, , drop = FALSE]); y <- e$pd_delta_other_kcal[complete]
  w <- as.numeric(weights(rep_eligible, type = "sampling")); wr <- as.matrix(weights(rep_eligible, type = "analysis"))
  full_response <- st_response_fit(R, complete, w)
  p <- full_response$probability
  probabilities <- data.frame(person_key = e$PERSON_KEY, complete = complete, probability = p, inverse_factor = 1 / p)
  diagnostics <- list(target_n = nrow(e), complete_n = sum(complete), target_key_sha256 = st_key_hash(e$PERSON_KEY),
    complete_key_sha256 = st_key_hash(e$PERSON_KEY[complete]), response_rank = full_response$rank,
    response_iterations = full_response$iterations, min_probability = min(p), max_probability = max(p),
    max_inverse_factor = max(1 / p), unweighted_response = mean(complete), weighted_response = weighted.mean(complete, w),
    parent_n = nrow(parent$variables), replicate_n = ncol(wr), df = survey::degf(rep_eligible),
    full_parent_df = survey::degf(rep_parent), mse = rep_eligible$mse, scale = rep_eligible$scale, rscales = rep_eligible$rscales,
    centering = "Full-sample point estimate", formula = deparse(response_formula), failed_replicates = 0L)
  # Positivity is evaluated for every eligible individual before adjusted outcome fitting.
  if (min(p) < cfg$minimum_probability || max(1 / p) > cfg$maximum_factor)
    return(list(status = "positivity_gate_failed_no_adjusted_outcome_fit", diagnostics = diagnostics, probabilities = probabilities))
  L <- st_L(colnames(X))
  full_beta <- list(baseline = st_wls(X, y, w[complete])$coefficients,
    response_weighted = st_wls(X, y, w[complete] / p[complete])$coefficients)
  full <- vapply(full_beta, function(b) sum(L * b), numeric(1))
  estimates <- matrix(NA_real_, ncol(wr), 2L, dimnames = list(NULL, names(full)))
  replicate_diagnostics <- vector("list", ncol(wr))
  for (j in seq_len(ncol(wr))) {
    result <- tryCatch({
      rw <- wr[, j]; rp <- st_response_fit(R, complete, rw)
      b <- st_wls(X, y, rw[complete])$coefficients
      a <- st_wls(X, y, rw[complete] / rp$probability[complete])$coefficients
      list(estimates = c(sum(L * b), sum(L * a)), response = rp)
    }, error = identity)
    failed <- inherits(result, "error")
    if (!failed) estimates[j, ] <- result$estimates
    replicate_diagnostics[[j]] <- data.frame(replicate = j, failed = failed,
      failure_reason = if (failed) conditionMessage(result) else "",
      response_rank = if (failed) NA_integer_ else result$response$rank,
      response_iterations = if (failed) NA_integer_ else result$response$iterations,
      minimum_probability = if (failed) NA_real_ else min(result$response$probability),
      maximum_factor = if (failed) NA_real_ else max(1 / result$response$probability))
    if (j %% 20L == 0L) progress(paste("JKn replicates", j, "of", ncol(wr)))
  }
  rd <- do.call(rbind, replicate_diagnostics); diagnostics$failed_replicates <- sum(rd$failed)
  weights_summary <- rbind(st_weight_diagnostics(w[complete], "original_complete_responders"),
    st_weight_diagnostics(w[complete] / p[complete], "response_adjusted_complete_responders"))
  out <- list(status = if (any(rd$failed)) "replicate_failure_no_inferential_output" else "complete",
    diagnostics = diagnostics, probabilities = probabilities, replicate_diagnostics = rd,
    replicate_estimates = estimates, internal_full_estimates = full,
    response_coefficients = full_response$coefficients, weights = weights_summary)
  if (any(rd$failed)) return(out)
  V <- survey::svrVar(estimates, scale = rep_eligible$scale, rscales = rep_eligible$rscales,
    mse = TRUE, coef = full)
  # Independently check package covariance against its explicit centered formula.
  centered <- sweep(estimates, 2, full, "-")
  explicit <- crossprod(centered, centered * rep_eligible$rscales) * rep_eligible$scale
  st_assert(st_near(V, explicit, 1e-10), "JKn scaling or MSE centering disagreement.")
  df <- survey::degf(rep_eligible); se <- sqrt(diag(V)); critical <- qt(.975, df)
  out$comparison <- data.frame(model = names(full), n = sum(complete), target_n = c(sum(complete), nrow(e)),
    eligible_population_n = nrow(e), target = c("complete_PHQ_pairs", "all_otherwise_eligible_pairs_under_conditional_response_exchangeability"),
    phq_ge10_n = sum(e$PHQ9_TOTAL[complete] >= 10), rank = ncol(X), df = df, estimate = as.numeric(full),
    SE = se, lower95 = full - critical * se, upper95 = full + critical * se, p_value = 2 * pt(-abs(full / se), df))
  out$covariance <- V
  out
}
