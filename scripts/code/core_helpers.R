# Exact current-analysis helper excerpts from the frozen source. See provenance/extractions.csv.
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
positive <- function(x) is.finite(x) & x > 0

read_project_csv <- function(path, nrows = -1L) {
  previous_ctype <- Sys.getlocale("LC_CTYPE")
  on.exit(suppressWarnings(Sys.setlocale("LC_CTYPE", previous_ctype)), add = TRUE)
  utf8_ctype <- suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
  assert(nzchar(utf8_ctype), "The native UTF-8 locale required for project CSV parsing is unavailable.")
  # Native read.csv(nrows=0) passes zero to scan's unlimited nmax setting.
  # Read one physical line ourselves when only the header is authorized.
  if (identical(as.integer(nrows), 0L)) {
    connection <- file(path, open = "rt", encoding = "UTF-8-BOM")
    header <- tryCatch(readLines(connection, n = 1L, warn = FALSE), finally = close(connection))
    assert(length(header) == 1L && nzchar(header), "Project CSV header is absent.")
    return(utils::read.csv(text = header, check.names = FALSE, na.strings = c("", "NA"), encoding = "UTF-8"))
  }
  utils::read.csv(path, check.names = FALSE, na.strings = c("", "NA"),
    fileEncoding = "UTF-8-BOM", encoding = "UTF-8", nrows = nrows)
}

parent_design <- function(d, weight = "WTMEC2YR") {
  keep <- positive(d[[weight]]) & is.finite(d$SDMVSTRA) & is.finite(d$SDMVPSU)
  z <- d[keep, , drop = FALSE]
  z$.weight <- z[[weight]] / 4
  z$.stratum <- z$SDMVSTRA
  z$.psu <- z$SDMVPSU
  assert(nrow(z) > 0L, "No positive-weight survey carrier remains.")
  survey::svydesign(ids = ~.psu, strata = ~.stratum, weights = ~.weight, data = z, nest = TRUE)
}

design_matrix <- function(fit, newdata) stats::model.matrix(stats::delete.response(stats::terms(fit)),
  newdata, contrasts.arg = fit$contrasts, xlev = fit$xlevels)

# Taylor influence for a population-standardized contrast includes variation in
# its weighted covariate distribution and its covariance with the fitted model.
# For the additive identity models the individual contrast is constant, so the
# empirical-distribution contribution cancels and L vcov(beta) L' is exact.
parent_score_covariance <- function(influence, domain_data, parent_data) {
  carrier <- matrix(0, nrow(parent_data), ncol(influence))
  carrier[match(domain_data$.row, parent_data$.row), ] <- influence
  psu <- interaction(parent_data$.stratum, parent_data$.psu, drop = TRUE)
  totals <- rowsum(carrier, psu, reorder = TRUE)
  ids <- unique(data.frame(key = as.character(psu), stratum = parent_data$.stratum))
  h <- ids$stratum[match(rownames(totals), ids$key)]
  covariance <- matrix(0, ncol(influence), ncol(influence))
  for (stratum in unique(h)) {
    values <- totals[h == stratum, , drop = FALSE]
    m <- nrow(values)
    assert(m > 1L, "The full parent design has a singleton PSU stratum.")
    centered <- sweep(values, 2, colMeans(values))
    covariance <- covariance + m / (m - 1) * crossprod(centered)
  }
  covariance
}

standardized_contrasts <- function(object, scenarios, combinations, return_influence = FALSE) {
  fit <- object$fit
  w <- object$data$.weight / sum(object$data$.weight)
  beta <- stats::coef(fit)
  V <- stats::vcov(fit)
  value <- numeric(length(scenarios))
  prediction <- matrix(NA_real_, nrow(object$data), length(scenarios))
  gradient <- matrix(NA_real_, length(scenarios), length(beta), dimnames = list(names(scenarios), names(beta)))
  for (i in seq_along(scenarios)) {
    X <- design_matrix(fit, scenarios[[i]])
    assert(identical(colnames(X), names(beta)), "Prediction matrix differs from the fitted coefficient space.")
    eta <- as.numeric(X %*% beta)
    mu <- fit$family$linkinv(eta)
    prediction[, i] <- mu
    value[i] <- sum(w * mu)
    gradient[i, ] <- colSums(X * (w * fit$family$mu.eta(eta)))
  }
  C <- as.matrix(combinations)
  assert(ncol(C) == length(scenarios), "Contrast scenarios do not align.")
  L <- C %*% gradient
  estimate <- as.numeric(C %*% value)
  names(estimate) <- rownames(C)
  conditional_covariance <- L %*% V %*% t(L)
  h <- prediction %*% t(C)
  empirical <- sweep(h, 2, estimate) * w
  X <- design_matrix(fit, object$data)
  eta <- as.numeric(X %*% beta)
  fitted <- fit$family$linkinv(eta)
  derivative <- fit$family$mu.eta(eta)
  variance <- fit$family$variance(fitted)
  bread <- crossprod(X, X * (object$data$.weight * derivative^2 / variance))
  score <- X * (object$data$.weight * as.numeric(stats::residuals(fit, type = "response")) * derivative / variance)
  influence <- empirical + score %*% solve(bread, t(L))
  covariance <- parent_score_covariance(influence, object$data, object$parent_design$variables)
  dimnames(covariance) <- list(names(estimate), names(estimate))
  out <- list(estimate = estimate, covariance = covariance, gradients = L,
    conditional_distribution_covariance = conditional_covariance,
    variance_method = "joint Taylor influence for fitted coefficients and weighted standardization distribution",
    df_complete = object$diagnostic$residual_df, diagnostic = object$diagnostic)
  if (return_influence) out$influence <- influence
  out
}

