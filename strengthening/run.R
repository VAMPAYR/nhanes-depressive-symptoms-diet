# Native Windows: Rscript --vanilla strengthening/run.R --preflight|--run|--verify [--output=PATH]
argv <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)), winslash = "/"))
root <- dirname(here)
args <- commandArgs(trailingOnly = TRUE)
mode <- intersect(args, c("--preflight", "--run", "--verify"))
if (length(mode) != 1L) stop("Specify exactly one of --preflight, --run, --verify.")
source(file.path(here, "config.R")); cfg <- strengthening_config
source(file.path(here, "analysis.R")); source(file.path(here, "tests.R"))
source(file.path(root, "scripts/config.R"))
if (nzchar(config$r_library)) .libPaths(c(file.path(root, "scripts", config$r_library), .libPaths()))
source(file.path(root, "scripts/public_validation.R"))
for (p in names(cfg$packages)) {
  st_assert(requireNamespace(p, quietly = TRUE), paste("Missing package", p))
  st_assert(packageVersion(p) == package_version(cfg$packages[[p]]), paste("Package version differs:", p))
}
st_assert(as.character(getRversion()) == cfg$expected_R, "R runtime version differs from frozen specification.")
for (file in c("paired_data.R", "master_builder.R", "paired_engine.R")) source(file.path(root, "scripts/code", file))
core <- new.env(parent = globalenv()); sys.source(file.path(root, "scripts/code/core_helpers.R"), envir = core)
data_root <- normalizePath(file.path(root, "scripts", config$data_root), winslash = "/")
plan <- file.path(here, "PLAN.md")
check_protected <- function() {
  public_check_manifest(root)
  public_compare(root, stages = c("original", "secondary"), write_report = FALSE)
  invisible(TRUE)
}
# Preflight validates inputs before original outputs exist on a clean checkout.
if (mode != "--preflight") check_protected()
outarg <- grep("^--output=", args, value = TRUE)
st_assert(length(outarg) <= 1L, "Use one output path.")
out <- if (length(outarg)) sub("^--output=", "", outarg) else file.path(here, cfg$output_relative)
out <- normalizePath(out, winslash = "/", mustWork = FALSE)
protected_roots <- normalizePath(file.path(root, "outputs", c("original", "secondary")), winslash = "/")
st_assert(!any(vapply(protected_roots, function(p) startsWith(tolower(paste0(out, "/")), tolower(paste0(p, "/"))), logical(1))), "Output path is protected.")
input_manifest <- read.csv(file.path(root, "scripts/raw_source_manifest.csv"), stringsAsFactors = FALSE)
input_manifest$path <- file.path(data_root, input_manifest$path)
extra <- do.call(rbind, lapply(names(pp_cycles()), function(cycle) data.frame(
  path = file.path(data_root, cycle, paste0("SMQ_", pp_cycles()[[cycle]], ".xpt")), sha256 = NA_character_, stringsAsFactors = FALSE)))
input_manifest <- rbind(input_manifest[c("path", "sha256")], extra)
input_manifest <- input_manifest[!duplicated(input_manifest$path), ]
st_assert(all(file.exists(input_manifest$path)), "A fixed input file is missing.")
actual <- vapply(input_manifest$path, st_hash, character(1))
known <- !is.na(input_manifest$sha256)
st_assert(all(actual[known] == input_manifest$sha256[known]), "Raw source differs from frozen original manifest.")
# SMQ inputs were audited and hashed before the analysis gate.
planning_manifest <- read.csv(file.path(here, "additional_raw_manifest.csv"), stringsAsFactors = FALSE)
smq <- grepl("/SMQ_", gsub("\\\\", "/", input_manifest$path))
mi <- match(normalizePath(input_manifest$path[smq], winslash = "/"), normalizePath(file.path(data_root, planning_manifest$path), winslash = "/"))
st_assert(!anyNA(mi) && all(actual[smq] == planning_manifest$sha256[mi]), "SMQ source differs from planning audit.")
input_manifest$sha256 <- actual
code_files <- c(list.files(here, pattern = "[.]R$", full.names = TRUE),
  file.path(root, "scripts/code", c("paired_data.R", "master_builder.R", "paired_engine.R", "core_helpers.R")),
  file.path(root, "scripts/config.R"), plan, file.path(here, "additional_raw_manifest.csv"))
code_manifest <- data.frame(path = normalizePath(code_files, winslash = "/"), sha256 = vapply(code_files, st_hash, character(1)))
if (mode == "--verify" || (mode == "--run" && file.exists(file.path(out, "COMPLETED.json")))) {
  receipt <- jsonlite::read_json(file.path(out, "COMPLETED.json"), simplifyVector = TRUE)
  st_assert(identical(receipt$status, "complete"), "Run lacks completion receipt.")
  saved_input <- read.csv(file.path(out, "input_manifest.csv"), stringsAsFactors = FALSE)
  saved_code <- read.csv(file.path(out, "code_manifest.csv"), stringsAsFactors = FALSE)
  st_assert(identical(saved_input$sha256, input_manifest$sha256) && identical(saved_code$sha256, code_manifest$sha256), "Run inputs or source have changed.")
  outputs <- read.csv(file.path(out, "output_manifest.csv"), stringsAsFactors = FALSE)
  st_assert(all(vapply(file.path(out, outputs$file), st_hash, character(1)) == outputs$sha256), "Output hash verification failed.")
  check_protected(); cat("Verified completed strengthening outputs and protected originals.\n"); quit(status = 0)
}
tests <- st_tests(core, cfg)
cat("Synthetic tests passed:", nrow(tests), "\n")
if (mode == "--preflight") { cat("Preflight passed. No participant outcome models fitted.\n"); quit(status = 0) }
st_assert(!dir.exists(out) || length(list.files(out, all.files = TRUE, no.. = TRUE)) == 0L,
  "Fresh --run requires an empty output directory; use --verify for an existing run.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(z, name) write.csv(z, file.path(out, name), row.names = FALSE, na = "")
write_json <- function(z, name) jsonlite::write_json(z, file.path(out, name), auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 16)
write_csv(input_manifest, "input_manifest.csv"); write_csv(code_manifest, "code_manifest.csv"); write_csv(tests, "synthetic_tests.csv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
write_json(list(status = "running", started_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE), distribution_spec_sha256 = st_hash(plan)), "STARTED.json")
run <- function() {
  cat("Rebuilding paired data from verified raw sources.\n")
  loaded <- st_load_data(data_root); d <- loaded$data
  counts <- c(parent = nrow(d), eligible = sum(d$st_eligible), complete = sum(d$st_eligible & d$st_complete),
    context = sum(d$st_eligible & d$st_complete & d$st_context_complete))
  st_assert(identical(as.integer(counts), as.integer(cfg$expected_n)), "Rebuilt feasibility counts disagree with frozen plan.")
  write_json(list(counts = as.list(counts), master_key_sha256 = st_key_hash(d$PERSON_KEY),
    primary_key_sha256 = st_key_hash(d$PERSON_KEY[d$st_eligible & d$st_complete]),
    context_key_sha256 = st_key_hash(d$PERSON_KEY[d$st_eligible & d$st_complete & d$st_context_complete]),
    input_policy = "Fresh raw-source rebuild; planning_master.rds is not an analysis input.",
    reproduction_scope = "Distribution rebuild; scientific functions unchanged from manuscript analysis."), "input_audit.json")
  saveRDS(loaded$builder_audit, file.path(out, "builder_audit.rds"))
  b <- st_missingness(d)
  for (name in names(b)) write_csv(b[[name]], paste0("B_", name, ".csv"))
  st_assert(sum(!d$st_dpq_record_present[d$st_eligible & !d$st_complete]) == 0 &&
    sum(d$st_all_items_NA[d$st_eligible & !d$st_complete]) == 100, "DPQ administration provenance differs from reviewed audit.")
  write_json(list(all_104_incomplete_have_source_record = TRUE, all_nine_NA_n = 100,
    administration = "DPQ documentation excludes proxy respondents; individual nonresponse reasons cannot be inferred from these records."), "B_administration_provenance.json")
  baseline <- list(); a <- list(); checks <- list()
  for (stage in c("primary", "quadratic")) {
    cat("Reconstructing frozen", stage, "model.\n")
    object <- pe_fit_object(d, core, stage)$object; f <- object$fit
    saved <- readRDS(file.path(root, "outputs/original", paste0("stage_", stage, ".rds")))$result
    # Original exported coefficient estimates preserve model-matrix ordering.
    estimate_column <- intersect(c("estimate", "Estimate"), names(saved$coefficients))
    st_assert(length(estimate_column) == 1L, "Saved coefficient estimate column is ambiguous.")
    coef_agree <- st_near(as.numeric(coef(f)), as.numeric(saved$coefficients[[estimate_column]]))
    covariance_agree <- st_near(vcov(f), saved$coefficient_covariance)
    st_assert(coef_agree && covariance_agree, paste("Frozen point/covariance reconstruction failed:", stage))
    row <- st_contrast(coef(f), vcov(f), f$df.residual, paste0(stage, "_original_Taylor"), nrow(object$data),
      sum(object$data$PHQ9_TOTAL >= 10), f$rank)
    baseline[[stage]] <- row
    checks[[stage]] <- data.frame(stage = stage, coefficient_max_absolute_difference = max(abs(coef(f) - saved$coefficients[[estimate_column]])),
      covariance_max_absolute_difference = max(abs(vcov(f) - saved$coefficient_covariance)), passed = TRUE)
    X <- model.matrix(f); a[[stage]] <- st_delete(X, object$data$pd_delta_other_kcal, object$data$.weight,
      object$data$SDMVSTRA, object$data$PERSON_KEY, row$SE, stage)
  }
  write_csv(do.call(rbind, checks), "baseline_reconstruction.csv")
  write_csv(do.call(rbind, baseline), "original_Taylor_comparisons.csv")
  for (name in c("cases", "strata", "summary")) write_csv(do.call(rbind, lapply(a, `[[`, name)), paste0("A_", name, ".csv"))
  cat("Fitting fixed same-sample context models.\n")
  context <- st_fit_context(d, core)
  write_csv(rbind(baseline$primary, context$comparison), "C_context_comparison.csv")
  saveRDS(context, file.path(out, "C_context_details.rds"))
  cat("Fitting fixed response model and full-parent JKn.\n")
  response <- st_response(d, core, cfg, function(message) cat(message, "\n"))
  write_json(response$diagnostics, "D_response_diagnostics.json")
  write_csv(response$probabilities, "D_response_probabilities.csv")
  if (!is.null(response$weights)) write_csv(response$weights, "D_weight_diagnostics.csv")
  if (!is.null(response$replicate_diagnostics)) write_csv(response$replicate_diagnostics, "D_replicate_diagnostics.csv")
  if (!is.null(response$replicate_estimates)) write_csv(data.frame(replicate = seq_len(nrow(response$replicate_estimates)), response$replicate_estimates), "D_replicate_estimates.csv")
  if (!is.null(response$comparison)) {
    st_assert(st_near(response$comparison$estimate[response$comparison$model == "baseline"], baseline$primary$estimate, 1e-8),
      "JKn baseline point differs from original Taylor point.")
    write_csv(response$comparison, "D_response_comparison.csv")
  }
  saveRDS(response, file.path(out, "D_response_details.rds"))
  check_protected()
  write_json(list(status = "complete", finished_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
    distribution_spec_sha256 = st_hash(plan), response_status = response$status,
    analysis_scope = "Fixed A-D; original six analyses remain unchanged; no C by D combination.",
    reference_aggregates_verified = TRUE), "COMPLETED.json")
  files <- list.files(out, full.names = TRUE)
  write_csv(data.frame(file = basename(files), sha256 = vapply(files, st_hash, character(1))), "output_manifest.csv")
  cat("Completed fixed strengthening analyses:", out, "\n")
}
tryCatch(run(), error = function(e) {
  write_json(list(status = "failed", error = conditionMessage(e), time_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE)), "FAILED.json")
  stop(conditionMessage(e), call. = FALSE)
})
