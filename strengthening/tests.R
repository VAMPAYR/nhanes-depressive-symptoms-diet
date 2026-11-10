# Synthetic tests only. No participant inputs or outcomes are accessed.
st_tests <- function(core, cfg) {
  set.seed(20261108)
  n <- 320L
  X <- cbind(1, pe_phq5 = rnorm(n), pe_phq5_squared = rnorm(n), nuisance = rnorm(n))
  y <- as.numeric(X %*% c(3, 2, -1, 4) + rnorm(n)); w <- exp(rnorm(n))
  z <- st_delete(X, y, w, rep(1:16, each = 20), seq_len(n), 1, "synthetic")
  L <- st_L(colnames(X))
  direct <- vapply(seq_len(n), function(i) sum(L * st_wls(X[-i, ], y[-i], w[-i])$coefficients), numeric(1))
  st_assert(st_near(z$cases$estimate, direct, 1e-9), "Synthetic deletion test failed.")
  d <- data.frame(RIDAGEYR = rep(40, 8), WTDR2D = 1, SDMVSTRA = 1, SDMVPSU = 1,
    pd_day1_valid = TRUE, pd_day2_valid = TRUE, pd_discordant = TRUE, pd_gap_days = 5,
    RIAGENDR = 1, RELEASE_ID = "2011-2012", pd_delta_other_kcal = 0, pd_order_sign = 1,
    PHQ9_TOTAL = c(0, 5, 10, NA, 0, 0, 0, 0), pd_base_eligible = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))
  for (field in paste0("pd_weekday", 2:7)) d[[field]] <- 0
  st_assert(all(st_eligible(d)), "Eligibility incorrectly depends on PHQ.")
  d$pd_day1_valid[1] <- FALSE; d$pd_gap_days[2] <- 2; d$pd_weekday2[3] <- 2
  st_assert(identical(st_eligible(d), c(FALSE, FALSE, FALSE, rep(TRUE, 5))), "Invalid pair eligibility accepted.")
  s <- data.frame(INDFMPIR = c(1, 2, 3, 4, 5, 6, NA, 0), DMDEDUC2 = c(1, 2, 3, 4, 5, 7, 9, NA),
    RIDRETH1 = c(1:5, 7, 9, NA), SMQ020 = c(2, 1, 1, 1, 7, 9, 1, 1), SMQ040 = c(NA, 3, 1, 2, 1, 3, 7, 9))
  s <- st_context(s)
  st_assert(identical(as.character(s$st_smoking), c("never", "former", "current", "current", rep(NA_character_, 4))), "Smoking special-code/skip test failed.")
  st_assert(identical(s$st_context_complete, c(rep(TRUE, 4), rep(FALSE, 4))), "Joint context completeness test failed.")
  # A full-parent design with a domain absent from selected parent rows.
  d <- data.frame(PERSON_KEY = seq_len(n), WTDR2D = w, SDMVSTRA = rep(1:32, each = 10),
    SDMVPSU = rep(rep(1:2, each = 5), 32), pe_phq5 = rnorm(n), pe_age10 = rnorm(n),
    pe_sex = factor(rep(1:2, length.out = n)), pe_cycle = factor(rep(1:4, length.out = n)),
    pd_order_sign = sample(c(-1, 1), n, TRUE), st_complete = runif(n) > .15,
    st_eligible = rep(c(TRUE, TRUE, TRUE, FALSE), length.out = n), pd_delta_other_kcal = rnorm(n) * 100)
  d$pe_age10_squared <- d$pe_age10^2; d$PHQ9_TOTAL <- pmax(0, d$pe_phq5 * 5)
  for (field in paste0("pd_weekday", 2:7)) d[[field]] <- sample(c(-1, 0, 1), n, TRUE)
  # Break deterministic cycle/sex aliasing while retaining all categories.
  d$pe_sex <- factor(sample(1:2, n, TRUE)); d$pe_cycle <- factor(sample(1:4, n, TRUE))
  response <- st_response(d, core, cfg)
  st_assert(identical(response$status, "complete") && response$diagnostics$failed_replicates == 0L,
    "Synthetic full-parent two-stage JKn failed.")
  st_assert(response$diagnostics$parent_n == n && response$diagnostics$target_n == sum(d$st_eligible),
    "Synthetic JKn did not preserve the full parent.")
  st_assert(all(is.finite(response$comparison$SE)) && all(response$comparison$SE > 0), "Synthetic JKn variance failed.")
  blocked <- cfg; blocked$minimum_probability <- 1
  stop_result <- st_response(d, core, blocked)
  st_assert(identical(stop_result$status, "positivity_gate_failed_no_adjusted_outcome_fit") &&
    is.null(stop_result$comparison) && is.null(stop_result$internal_full_estimates), "Positivity gate ran outcome fitting.")
  # Same-sample baseline and extended model, with deliberate context missingness.
  d$INDFMPIR <- runif(n, 0, 5); d$DMDEDUC2 <- sample(1:5, n, TRUE); d$RIDRETH1 <- sample(1:5, n, TRUE)
  d$SMQ020 <- sample(1:2, n, TRUE); d$SMQ040 <- sample(1:3, n, TRUE); d$INDFMPIR[1:12] <- NA
  dc <- st_context(d); context <- st_fit_context(dc, core)
  expected <- sum(dc$st_eligible & dc$st_complete & dc$st_context_complete)
  st_assert(all(context$comparison$n == expected), "Same-sample context test failed.")
  data.frame(test = c("all_exact_deletions_against_direct_WLS", "PHQ_independent_eligibility", "smoking_codes_and_skips",
    "joint_context_completeness", "full_parent_two_stage_JKn", "positivity_gate_precedes_outcomes", "same_sample_context_models"), passed = TRUE)
}
