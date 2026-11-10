#!/usr/bin/env Rscript
# Separate optional entry point; the original three-fit entry point is unchanged.
entry <- grep("^--file=", commandArgs(), value = TRUE)
if (length(entry) != 1L) stop("Run this file with Rscript.", call. = FALSE)
package_root <- normalizePath(dirname(sub("^--file=", "", entry)), winslash = "/", mustWork = TRUE)
source(file.path(package_root, "config.R"), local = TRUE)
resolve <- function(path) {
  if (!grepl("^([A-Za-z]:|/|\\\\)", path)) path <- file.path(package_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}
if (nzchar(config$r_library)) .libPaths(c(resolve(config$r_library), .libPaths()))
for (name in c("core_helpers.R", "checkpoint_helpers.R", "paired_data.R", "paired_engine.R",
  "master_builder.R", "secondary_engine.R")) source(file.path(package_root, "code", name), local = TRUE)
options(stringsAsFactors = FALSE, survey.lonely.psu = "fail", warn = 1)
data_root <- resolve(config$data_root)
output_root <- resolve(if (is.null(config$secondary_output_root)) "secondary_output" else config$secondary_output_root)
assert(!identical(output_root, resolve(config$output_root)), "Secondary output must differ from primary output.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) args <- "--preflight"
assert(length(args) == 1L && args %in% c("--preflight", "--run", "--verify"), "Use --preflight, --run, or --verify.")

secondary_preflight <- function() {
  assert(as.character(getRversion()) == config$expected_R,
    paste("Use R", config$expected_R))
  for (name in names(config$packages)) assert(requireNamespace(name, quietly = TRUE) &&
    utils::compareVersion(as.character(utils::packageVersion(name)), config$packages[[name]]) == 0,
    paste("Missing or different installed package:", name))
  manifests <- lapply(c("active_code_manifest.csv", "secondary_code_manifest.csv", "raw_source_manifest.csv"),
    function(name) utils::read.csv(file.path(package_root, name), stringsAsFactors = FALSE))
  for (i in seq_along(manifests)) {
    manifest <- manifests[[i]]
    assert(!anyDuplicated(manifest$path) &&
      all(!grepl("(^[/\\\\]|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))", manifest$path)), "Unsafe manifest path.")
    root <- if (i == 3L) data_root else package_root
    for (j in seq_len(nrow(manifest))) assert(file.exists(file.path(root, manifest$path[j])) &&
      identical(oa_sha(file.path(root, manifest$path[j])), manifest$sha256[j]),
      paste("Source missing or changed:", manifest$path[j]))
  }
  assert(identical(pd_stage_specs()$stage_id, c("primary", "all_intervals", "quadratic")), "Core fit set changed.")
  message("Secondary preflight passed; no participant records parsed; exactly three optional outcome fits specified.")
  list(package = "paired_current_20261102_secondary", R = as.character(getRversion()), packages = config$packages,
    core_code = manifests[[1]], secondary_code = manifests[[2]], raw = manifests[[3]][c("path", "sha256")],
    stages = se_stages(), prior_primary_results_seen = TRUE, imputation = FALSE, model_selection = FALSE,
    inference = "Stored svyglm residual df; Holm across two interaction tests")
}
identity <- secondary_preflight()
if (args == "--preflight") quit(status = 0L)
marker_path <- file.path(output_root, "COMPLETED.json")
secondary_verify <- function() {
  marker <- jsonlite::read_json(marker_path, simplifyVector = TRUE)
  assert(isTRUE(marker$completed) && identical(marker$identity_sha256, digest::digest(identity, algo = "sha256")),
    "Secondary completion identity changed.")
  for (name in names(marker$outputs)) assert(identical(oa_sha(file.path(output_root, name)), marker$outputs[[name]]),
    paste("Secondary output changed:", name))
  message("Secondary completion identity and all output hashes verified.")
  invisible(TRUE)
}
if (args == "--verify") { secondary_verify(); quit(status = 0L) }

secondary_run <- function() {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(marker_path)) { secondary_verify(); message("No secondary stages repeated."); return(invisible(TRUE)) }
  lease <- oa_acquire_lock(output_root)
  on.exit(oa_release_lock(lease), add = TRUE)
  started <- proc.time()[["elapsed"]]
  status <- function(stage, detail = "") oa_json(list(stage = stage, detail = detail,
    specified_outcome_fits = 3L, primary_outcome_fits = 0L, elapsed_seconds = proc.time()[["elapsed"]] - started),
    file.path(output_root, "STATUS.json"))
  bounded <- function(compute) {
    remaining <- 600 - (proc.time()[["elapsed"]] - started)
    assert(remaining > 0, "Secondary ten-minute runtime guard expired.")
    setTimeLimit(cpu = Inf, elapsed = min(config$stage_elapsed_limit_seconds, remaining), transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
    compute()
  }
  checkpoint <- function(id, compute) {
    status(id)
    result <- oa_stage_checkpoint(file.path(output_root, paste0("stage_", id, ".rds")),
      list(identity = identity, stage = id), function() bounded(compute))
    message(if (result$reused) "Reused: " else "Saved: ", id)
    result$result
  }
  tryCatch({
    binding_path <- file.path(output_root, "run_identity.rds")
    if (file.exists(binding_path)) assert(identical(readRDS(binding_path), identity), "Secondary checkpoint identity changed.")
    else oa_save_once(identity, binding_path)
    # Save the already fixed specification before participant parsing or new fits.
    spec_path <- file.path(output_root, "FIXED_SPEC_BEFORE_EXECUTION.txt")
    spec <- readLines(file.path(package_root, "SECONDARY_SPEC.md"), warn = FALSE)
    if (file.exists(spec_path)) assert(identical(readLines(spec_path, warn = FALSE), spec), "Saved secondary specification changed.")
    else writeLines(spec, spec_path)
    status("raw_rebuild")
    master <- bounded(function() mb_build(data_root))
    built <- bounded(function() pp_build_data(master$data, data_root))
    prepared <- se_prepare(built$data)
    rm(master, built)
    cohort <- list(primary_n = sum(prepared$data$.se_domain), ordered_person_key_sha256 = prepared$key_sha256,
      hash_method = "digest SHA256 of original ordered character PERSON_KEY vector; keys are not saved",
      masks = prepared$counts)
    cohort_path <- file.path(output_root, "cohort_identity.rds")
    if (file.exists(cohort_path)) assert(identical(readRDS(cohort_path), cohort), "Secondary cohort identity changed.")
    else oa_save_once(cohort, cohort_path)
    fits <- list()
    for (id in se_stages()) {
      stage <- id
      fits[[stage]] <- checkpoint(stage, function() se_fit(prepared$data, stage))
    }
    patterns <- checkpoint("recall_pattern_descriptions", function() se_patterns(prepared$data))
    tests <- do.call(rbind, lapply(fits[c("PHQ_order", "PHQ_sex")], `[[`, "test"))
    tests$p_Holm_two_interactions <- stats::p.adjust(tests$p_value, method = "holm")
    tagged_covariance <- function(result, stage) {
      V <- result$contrast_covariance
      out <- expand.grid(estimand_row = rownames(V), estimand_column = colnames(V), stringsAsFactors = FALSE)
      out$covariance <- as.vector(V); out$stage_id <- stage; out
    }
    tables <- list(secondary_contrasts = do.call(rbind, lapply(fits, `[[`, "table")),
      secondary_coefficients = do.call(rbind, lapply(fits, `[[`, "coefficients")),
      interaction_tests = tests, recall_patterns_by_sex_PHQ = patterns$table,
      subgroup_support = oa_bind_tables(lapply(fits, `[[`, "support")),
      contrast_covariance = do.call(rbind, lapply(names(fits), function(id) tagged_covariance(fits[[id]], id))),
      analysis_domain_flow = prepared$counts)
    oa_publish_tables(tables, output_root)
    details <- list(cohort = cohort, diagnostics = lapply(fits, `[[`, "diagnostic"),
      functional_checks = lapply(fits, `[[`, "validation"), descriptive_metadata = patterns$metadata,
      prior_primary_results_seen = TRUE, secondary_only_outcome_fits = se_stages(),
      primary_runner_unchanged = TRUE, no_participant_records_saved = TRUE)
    details_path <- file.path(output_root, "analysis_details.json")
    if (file.exists(details_path)) assert(identical(jsonlite::toJSON(details, auto_unbox = TRUE, na = "null", null = "null", digits = NA),
      jsonlite::toJSON(jsonlite::read_json(details_path, simplifyVector = FALSE), auto_unbox = TRUE, na = "null", null = "null", digits = NA)),
      "Saved secondary details changed.") else oa_json(details, details_path)
    assert(identical(secondary_preflight(), identity), "Secondary inputs changed during execution.")
    files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
    files <- files[!dir.exists(files) & !grepl("/(RUNNING_LOCK/|STATUS\\.json$)", files)]
    relative <- substring(files, nchar(output_root) + 2L)
    oa_json(list(completed = TRUE, identity_sha256 = digest::digest(identity, algo = "sha256"),
      outcome_fits = se_stages(), primary_fits_executed = 0L,
      elapsed_seconds = proc.time()[["elapsed"]] - started,
      outputs = setNames(as.list(vapply(files, oa_sha, "")), relative)), marker_path)
    status("completed")
    message("Completed exactly three secondary fits and unadjusted recall-pattern descriptions.")
    invisible(TRUE)
  }, error = function(e) { status("failed", conditionMessage(e)); stop(conditionMessage(e), call. = FALSE) },
  interrupt = function(e) { status("interrupted", "Checkpoints preserved; no automatic retry."); stop("Interrupted.", call. = FALSE) })
}
secondary_run()
