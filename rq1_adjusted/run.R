# Native Windows: Rscript --vanilla rq1_adjusted/run.R --preflight|--run|--verify [--output=PATH]
# Fourth stage: order- and weekday-adjusted RQ1 energy partition (frozen SPEC.md). Post hoc sensitivity analysis
# specified after the original results were known and before fitting; the original, secondary and strengthening
# outputs are protected and never rewritten.
argv <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)), winslash = "/"))
root <- dirname(here)
args <- commandArgs(trailingOnly = TRUE)
rq_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
mode <- intersect(args, c("--preflight", "--run", "--verify"))
rq_assert(length(mode) == 1L, "Specify exactly one of --preflight, --run, --verify.")

rq_config <- list(
  spec_sha256 = "9293e87e636aca197817af1a57a34318083c7e490da8273238d71050b851a448",
  expected_R = "4.5.1", packages = c(survey = "4.5", foreign = "0.8.90", digest = "0.6.37", jsonlite = "2.0.0"),
  output_relative = "../outputs/rq1_adjusted/results", expected_n = 2320L,
  reference_tolerance = 1e-8, table2_tolerance = 1e-6)
cfg <- rq_config
for (p in names(cfg$packages)) {
  rq_assert(requireNamespace(p, quietly = TRUE), paste("Missing package", p))
  rq_assert(packageVersion(p) == package_version(cfg$packages[[p]]), paste("Package version differs:", p))
}
rq_assert(as.character(getRversion()) == cfg$expected_R, "R runtime version differs from the frozen specification.")
rq_hash <- function(path) digest::digest(file = path, algo = "sha256")

source(file.path(root, "scripts/config.R"))
if (nzchar(config$r_library)) .libPaths(c(file.path(root, "scripts", config$r_library), .libPaths()))
source(file.path(root, "scripts/public_validation.R"))
for (file in c("paired_data.R", "master_builder.R", "paired_engine.R")) source(file.path(root, "scripts/code", file))
core <- new.env(parent = globalenv()); sys.source(file.path(root, "scripts/code/core_helpers.R"), envir = core)
data_root <- normalizePath(file.path(root, "scripts", config$data_root), winslash = "/")
spec <- file.path(here, "SPEC.md")
rq_assert(identical(rq_hash(spec), cfg$spec_sha256), "Frozen specification hash differs.")

check_protected <- function() {
  public_check_manifest(root)
  public_compare(root, stages = c("original", "secondary"), write_report = FALSE)
  invisible(TRUE)
}
# Preflight validates inputs before original outputs exist on a clean checkout.
if (mode != "--preflight") check_protected()

outarg <- grep("^--output=", args, value = TRUE)
rq_assert(length(outarg) <= 1L, "Use one output path.")
out <- if (length(outarg)) sub("^--output=", "", outarg) else file.path(here, cfg$output_relative)
out <- normalizePath(out, winslash = "/", mustWork = FALSE)
protected_roots <- normalizePath(file.path(root, "outputs", c("original", "secondary", "strengthening")), winslash = "/")
rq_assert(!any(vapply(protected_roots, function(p) startsWith(tolower(paste0(out, "/")), tolower(paste0(p, "/"))), logical(1))),
  "Output path is protected.")

input_manifest <- read.csv(file.path(root, "scripts/raw_source_manifest.csv"), stringsAsFactors = FALSE)[c("path", "sha256")]
input_manifest$path <- file.path(data_root, input_manifest$path)
rq_assert(all(file.exists(input_manifest$path)), "A fixed raw input file is missing.")
actual <- vapply(input_manifest$path, rq_hash, character(1))
rq_assert(all(actual == input_manifest$sha256), "Raw source differs from the frozen original manifest.")
code_files <- c(file.path(here, c("run.R", "SPEC.md")),
  file.path(root, "reference/rq1_adjusted", c("rq1_adjusted_results.csv", "rq1_adjusted_summary.json")),
  file.path(root, "scripts/code", c("paired_data.R", "master_builder.R", "paired_engine.R", "core_helpers.R")),
  file.path(root, "scripts/config.R"))
code_manifest <- data.frame(path = normalizePath(code_files, winslash = "/"), sha256 = vapply(code_files, rq_hash, character(1)))

components <- c(beverage_ethanol_kcal = "Ethanol energy in coded beverages",
  beverage_nonethanol_kcal = "Calculated coded-beverage energy remainder",
  linked_addition_kcal = "Energy in linked beverage additions",
  other_kcal = "Energy from other reported foods and beverages", total_kcal = "Total reported energy")
rq_formula <- function(response) stats::reformulate(c("pd_order_sign", paste0("pd_weekday", 2:7)), response = response)
rq_fit <- function(design, response) {
  fit <- survey::svyglm(rq_formula(response), design = design, family = stats::gaussian(link = "identity"),
    control = stats::glm.control(epsilon = 1e-10, maxit = 100))
  rq_assert(isTRUE(fit$converged) && all(is.finite(stats::coef(fit))), paste("Adjusted fit failed:", response))
  fit
}

# Synthetic checks: exact recovery of a known intercept, additivity of component intercepts.
rq_tests <- function() {
  set.seed(20261120)
  n <- 400L
  d <- data.frame(SDMVSTRA = rep(1:20, each = 20), SDMVPSU = rep(rep(1:2, each = 10), 20), WTDR2D = stats::runif(n, 1, 3))
  d$pd_order_sign <- sample(c(-1, 1), n, replace = TRUE)
  drink <- sample(1:7, n, replace = TRUE); comparison <- sample(1:7, n, replace = TRUE)
  for (k in 2:7) d[[paste0("pd_weekday", k)]] <- as.numeric(drink == k) - as.numeric(comparison == k)
  d$pd_delta_a <- 100 + 5 * d$pd_order_sign
  d$pd_delta_b <- 50 - 40 * d$pd_order_sign + 10 * d$pd_weekday7 + stats::rnorm(n)
  d$pd_delta_t <- d$pd_delta_a + d$pd_delta_b
  design <- core$parent_design(d, weight = "WTDR2D")
  b <- vapply(c("pd_delta_a", "pd_delta_b", "pd_delta_t"), function(v) stats::coef(rq_fit(design, v))[["(Intercept)"]], numeric(1))
  tests <- data.frame(test = c("exact known intercept recovered", "noisy intercept within 1 kcal of truth",
    "component intercepts sum to the total intercept"),
    passed = c(abs(b[[1]] - 100) < 1e-8, abs(b[[2]] - 50) < 1, abs(b[[1]] + b[[2]] - b[[3]]) < 1e-8))
  rq_assert(all(tests$passed), "A synthetic check failed.")
  tests
}

if (mode == "--verify" || (mode == "--run" && file.exists(file.path(out, "COMPLETED.json")))) {
  receipt <- jsonlite::read_json(file.path(out, "COMPLETED.json"), simplifyVector = TRUE)
  rq_assert(identical(receipt$status, "complete"), "Run lacks a completion receipt.")
  saved_input <- read.csv(file.path(out, "input_manifest.csv"), stringsAsFactors = FALSE)
  saved_code <- read.csv(file.path(out, "code_manifest.csv"), stringsAsFactors = FALSE)
  rq_assert(identical(saved_input$sha256, unname(actual)) && identical(saved_code$sha256, unname(code_manifest$sha256)),
    "Run inputs or source have changed.")
  outputs <- read.csv(file.path(out, "output_manifest.csv"), stringsAsFactors = FALSE)
  rq_assert(all(vapply(file.path(out, outputs$file), rq_hash, character(1)) == outputs$sha256), "Output hash verification failed.")
  check_protected(); public_compare(root, stages = "rq1_adjusted", write_report = FALSE)
  cat("Verified completed RQ1 adjusted outputs and fixed references.\n"); quit(status = 0)
}
tests <- rq_tests()
cat("Synthetic tests passed:", nrow(tests), "\n")
if (mode == "--preflight") { cat("Preflight passed. No participant outcome models fitted.\n"); quit(status = 0) }

rq_assert(!dir.exists(out) || length(list.files(out, all.files = TRUE, no.. = TRUE)) == 0L,
  "A fresh --run requires an empty output directory; use --verify for an existing run.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(z, name) write.csv(z, file.path(out, name), row.names = FALSE, na = "")
write_json <- function(z, name) jsonlite::write_json(z, file.path(out, name), auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 16)
input_manifest$sha256 <- unname(actual)
write_csv(input_manifest, "input_manifest.csv"); write_csv(code_manifest, "code_manifest.csv"); write_csv(tests, "synthetic_tests.csv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
write_json(list(status = "running", started_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), spec_sha256 = cfg$spec_sha256), "STARTED.json")

cat("Rebuilding paired data from verified raw sources.\n")
paired <- pp_build_data(mb_build(data_root)$data, data_root)
prepared <- pe_prepare(paired$data)
d <- prepared$data
d$.rq_domain <- prepared$masks$primary_3_10
rq_assert(sum(d$.rq_domain) == cfg$expected_n, "Primary domain size differs from 2,320.")
parent <- core$parent_design(d, weight = "WTDR2D")
domain <- subset(parent, .rq_domain)
design_df <- survey::degf(parent)

# Gate: reproduce Table 2 (unadjusted oriented differences) before reading any adjusted model.
saved <- read.csv(file.path(root, "scripts/reference_outputs/paired/weighted_energy_decomposition.csv"), stringsAsFactors = FALSE)
saved <- saved[saved$window == "primary_days3to10" & saved$PHQ_band == "All PHQ scores" & saved$state == "oriented_difference", ]
gate <- do.call(rbind, lapply(names(components), function(component) {
  m <- survey::svymean(stats::as.formula(paste0("~pd_delta_", component)), domain)
  est <- as.numeric(stats::coef(m)); se <- as.numeric(survey::SE(m)); q <- stats::qt(.975, design_df)
  ref <- saved[saved$variable == paste0(".pp_", component), ]
  rq_assert(nrow(ref) == 1L, paste("Table 2 reference row missing:", component))
  data.frame(component = components[[component]], estimate = est, lower = est - q * se, upper = est + q * se,
    table2_estimate = ref$estimate, table2_lower = ref$ci_lower, table2_upper = ref$ci_upper,
    max_absolute_difference = max(abs(c(est, est - q * se, est + q * se) - c(ref$estimate, ref$ci_lower, ref$ci_upper))))
}))
gate$passed <- gate$max_absolute_difference < cfg$table2_tolerance
write_csv(gate, "gate_table2_reproduction.csv")
rq_assert(all(gate$passed), "Gate failed: Table 2 not reproduced. No adjusted model is read.")

cat("Fitting the five fixed adjusted models.\n")
results <- do.call(rbind, lapply(names(components), function(component) {
  fit <- rq_fit(domain, paste0("pd_delta_", component))
  b <- stats::coef(fit); V <- stats::vcov(fit); q <- stats::qt(.975, fit$df.residual)
  data.frame(component = components[[component]], adjusted = b[["(Intercept)"]], se = sqrt(V["(Intercept)", "(Intercept)"]),
    lower = b[["(Intercept)"]] - q * sqrt(V["(Intercept)", "(Intercept)"]),
    upper = b[["(Intercept)"]] + q * sqrt(V["(Intercept)", "(Intercept)"]),
    order_coef = b[["pd_order_sign"]],
    order_lower = b[["pd_order_sign"]] - q * sqrt(V["pd_order_sign", "pd_order_sign"]),
    order_upper = b[["pd_order_sign"]] + q * sqrt(V["pd_order_sign", "pd_order_sign"]),
    residual_df = fit$df.residual, n = nrow(fit$data), row.names = NULL)
}))
identity_gap <- sum(results$adjusted[1:4]) - results$adjusted[5]
rq_assert(abs(identity_gap) < 1e-6, "Adjusted components do not sum to the adjusted total.")
write_csv(results, "rq1_adjusted_results.csv")

reference <- read.csv(file.path(root, "reference/rq1_adjusted/rq1_adjusted_results.csv"), stringsAsFactors = FALSE)
fields <- c("adjusted", "se", "lower", "upper", "order_coef", "order_lower", "order_upper", "residual_df", "n")
agreement <- data.frame(field = fields, max_absolute_difference = vapply(fields, function(f)
  max(abs(as.numeric(results[[f]]) - as.numeric(reference[[f]]))), numeric(1)))
agreement$passed <- agreement$max_absolute_difference < cfg$reference_tolerance
write_csv(agreement, "reference_agreement.csv")
rq_assert(all(agreement$passed), "Raw-source rebuild disagrees with the first computation.")

beverage <- sum(results$adjusted[1:3])
write_json(list(n = cfg$expected_n, design_df = design_df, residual_df = results$residual_df[1],
  identity_gap_kcal = identity_gap, beverage_sum_kcal = beverage, beverage_share = beverage / results$adjusted[5],
  input_policy = "Fresh raw-source rebuild; planning_master.rds is not an input.",
  status_note = "Post hoc sensitivity analysis, specified before fitting after the original results were known."),
  "rq1_adjusted_summary.json")
check_protected()
public_compare(root, stages = "rq1_adjusted", write_report = FALSE)
write_json(list(status = "complete", finished_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), spec_sha256 = cfg$spec_sha256,
  analysis_scope = "Five fixed adjusted component models; Table 2 and the six existing models unchanged.",
  reference_aggregates_verified = TRUE), "COMPLETED.json")
files <- setdiff(list.files(out), "output_manifest.csv")
write_csv(data.frame(file = files, sha256 = vapply(file.path(out, files), rq_hash, character(1))), "output_manifest.csv")
cat("RQ1 adjusted stage complete:", out, "\n")
