#!/usr/bin/env Rscript
# One entry point. No participant data are read by --preflight.
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
  "master_builder.R", "manuscript_descriptives.R", "reporting.R")) source(file.path(package_root, "code", name), local = TRUE)
options(stringsAsFactors = FALSE, survey.lonely.psu = "fail", warn = 1)
data_root <- resolve(config$data_root)
output_root <- resolve(config$output_root)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) args <- "--preflight"
assert(length(args) == 1L && args %in% c("--preflight", "--run", "--verify"), "Use --preflight, --run, or --verify.")

portable_preflight <- function() {
  assert(as.character(getRversion()) == config$expected_R,
    paste("Use R", config$expected_R))
  for (name in names(config$packages)) assert(requireNamespace(name, quietly = TRUE) &&
    utils::compareVersion(as.character(utils::packageVersion(name)), config$packages[[name]]) == 0,
    paste("Missing or different installed package:", name, config$packages[[name]]))
  raw <- utils::read.csv(file.path(package_root, "raw_source_manifest.csv"), stringsAsFactors = FALSE)
  code <- utils::read.csv(file.path(package_root, "active_code_manifest.csv"), stringsAsFactors = FALSE)
  for (manifest in list(raw, code)) assert(!anyDuplicated(manifest$path) &&
    all(!grepl("(^[/\\\\]|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))", manifest$path)), "Unsafe manifest path.")
  for (i in seq_len(nrow(code))) assert(identical(pm_sha(file.path(package_root, code$path[i])), code$sha256[i]),
    paste("Package code changed:", code$path[i]))
  for (i in seq_len(nrow(raw))) {
    path <- file.path(data_root, raw$path[i])
    assert(file.exists(path) && identical(pm_sha(path), raw$sha256[i]), paste("Raw source missing or changed:", path))
  }
  assert(is.numeric(config$stage_elapsed_limit_seconds) && config$stage_elapsed_limit_seconds > 0,
    "Invalid cooperative stage time limit.")
  identity <- list(package = "paired_current_20261102", R = as.character(getRversion()), packages = as.list(config$packages),
    code = code, raw_sources = raw[c("path", "sha256")], fits = c("primary", "all_intervals", "quadratic"),
    inference = "survey_default", imputation = FALSE, model_selection = FALSE)
  message("Preflight passed: ", nrow(raw), " raw sources (", round(sum(raw$bytes) / 1024^2, 1),
    " MiB); fixed three-fit workload; no participant records read.")
  identity
}

portable_verify <- function(marker) {
  assert(isTRUE(marker$completed), "Output completion marker is not complete.")
  assert(identical(marker$identity_sha256, digest::digest(identity, algo = "sha256")), "Completed input/code identity differs.")
  for (name in names(marker$outputs)) assert(identical(pm_sha(file.path(output_root, name)), marker$outputs[[name]]),
    paste("Completed output changed:", name))
  message("Completed outputs and their source identity verified.")
  invisible(TRUE)
}

compare_reference <- function() {
  results <- list()
  for (section in c("paired", "descriptives")) {
    baseline <- file.path(package_root, "reference_outputs", section)
    fresh <- if (section == "paired") output_root else file.path(output_root, "descriptives")
    for (file in list.files(baseline, pattern = "\\.csv$")) {
      old <- utils::read.csv(file.path(baseline, file), stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
      now <- utils::read.csv(file.path(fresh, file), stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
      matched <- isTRUE(all.equal(now, old, tolerance = 1e-9, check.attributes = FALSE))
      results[[length(results) + 1L]] <- data.frame(section = section, file = file,
        rows = nrow(now), columns = ncol(now), numeric_tolerance = 1e-9, matched = matched)
      assert(matched, paste("Reference aggregate differs:", section, file))
    }
  }
  figure <- file.path(output_root, "figures/paired_phq_profiles.png")
  figure_match <- identical(pm_sha(figure), pm_sha(file.path(package_root, "reference_outputs/figures/paired_phq_profiles.png")))
  list(tables = do.call(rbind, results), figure_bytes_identical = figure_match,
    figure_scope = "Same saved point estimates/intervals and drawing specification; byte equality also depends on graphics device and font installation")
}

identity <- portable_preflight()
if (args == "--preflight") quit(status = 0L)
marker_path <- file.path(output_root, "COMPLETED.json")
if (args == "--verify") {
  portable_verify(jsonlite::read_json(marker_path, simplifyVector = TRUE))
  comparison <- compare_reference()
  print(comparison$tables)
  message("Publication figure byte match: ", comparison$figure_bytes_identical)
  quit(status = 0L)
}

run_portable <- function() {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(marker_path)) {
    portable_verify(jsonlite::read_json(marker_path, simplifyVector = TRUE))
    message("The fixed analysis is complete. No stages repeated.")
    return(invisible(TRUE))
  }
  lease <- oa_acquire_lock(output_root)
  on.exit(oa_release_lock(lease), add = TRUE)
  identity_path <- file.path(output_root, "run_identity.rds")
  if (file.exists(identity_path)) assert(identical(readRDS(identity_path), identity), "Checkpoint source identity differs.")
  else oa_save_once(identity, identity_path)
  started <- proc.time()[["elapsed"]]
  status <- function(stage, detail = "") oa_json(list(stage = stage, detail = detail,
    elapsed_seconds = proc.time()[["elapsed"]] - started, specified_outcome_fits = 3L,
    imputation = FALSE, model_selection = FALSE), file.path(output_root, "STATUS.json"))
  bounded <- function(compute) {
    # The total guard includes all stages and stops at R's next interrupt point.
    remaining <- 600 - (proc.time()[["elapsed"]] - started)
    assert(remaining > 0, "The ten-minute total cooperative runtime guard expired.")
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
    status("master_and_paired_construction")
    master <- bounded(function() mb_build(data_root))
    utils::write.csv(master$join_audit, file.path(output_root, "master_join_audit.csv"), row.names = FALSE, na = "")
    oa_json(master$metadata, file.path(output_root, "master_build_metadata.json"))
    built <- bounded(function() pp_build_data(master$data, data_root))
    rm(master)
    domains <- pd_domains(built$data)
    assert(sum(domains$masks$primary_3_10) == 2320L && sum(domains$masks$all_intervals) == 2701L,
      "Canonical model cohorts changed.")
    specs <- pd_stage_specs()
    assert(identical(specs$stage_id, c("primary", "all_intervals", "quadratic")), "The fixed fit set changed.")
    core <- environment(parent_design)
    paired <- environment(pd_fit_stage)
    fits <- list()
    for (id in specs$stage_id) {
      stage_id <- id
      fits[[id]] <- checkpoint(id, function() {
        result <- pd_fit_stage(built$data, core, stage_id, config = list(residual_df_policy = "survey_default"))
        result$table <- result$profiles
        result
      })
    }
    context <- checkpoint("paired_descriptions", function() {
      result <- pp_descriptives(built$data, core)
      result$source_audit <- built$audit
      result$table <- result$flow
      result
    })
    tagged <- function(key) oa_bind_tables(lapply(names(fits), function(id) {
      z <- fits[[id]][[key]]; z$stage_id <- id; z
    }))
    tables <- list(PHQ_profiles = tagged("profiles"), PHQ_contrasts = tagged("contrasts"),
      model_coefficients = tagged("coefficients"), analysis_domain_flow = domains$counts,
      primary_estimate = fits$primary$primary_estimate,
      primary_test = data.frame(fits$primary$primary_wald, stringsAsFactors = FALSE))
    for (name in names(context)) if (is.data.frame(context[[name]]) && name != "table") tables[[name]] <- context[[name]]
    tables$model_diagnostics <- oa_bind_tables(lapply(names(fits), function(id) {
      z <- fits[[id]]$diagnostic
      scalars <- z[vapply(z, function(x) length(x) == 1 && is.atomic(x), logical(1))]
      data.frame(stage_id = id, scalars, check.names = FALSE, stringsAsFactors = FALSE)
    }))
    tables$residual_summaries <- oa_bind_tables(lapply(names(fits), function(id) {
      z <- fits[[id]]$residual_diagnostics$grouped; z$stage_id <- id; z
    }))
    tables$leverage_summary <- oa_bind_tables(lapply(names(fits), function(id) {
      z <- fits[[id]]$residual_diagnostics$leverage; z$stage_id <- id; z
    }))
    for (level in c("summary", "strata", "top_psus")) tables[[paste0("influence_", level)]] <-
      oa_bind_tables(lapply(names(fits), function(id) {
        z <- fits[[id]]$influence_contributions[[level]]; z$stage_id <- id; z
      }))
    status("aggregate_exports")
    oa_publish_tables(tables, output_root)
    oa_json(list(question = "PHQ heterogeneity in paired other-food/beverage energy differences",
      fixed_stages = specs, primary_test = fits$primary$primary_wald,
      design_diagnostics = lapply(fits, `[[`, "diagnostic"), source_audit = built$audit,
      residual_diagnostics = lapply(fits, `[[`, "residual_diagnostics"), context_metadata = context$metadata),
      file.path(output_root, "analysis_details.json"))
    descriptive_dir <- file.path(output_root, "descriptives")
    dir.create(descriptive_dir, showWarnings = FALSE)
    descriptive <- checkpoint("manuscript_descriptions", function() {
      result <- list(tables = rp_descriptives(built$data, core, paired))
      result$table <- result$tables$table1_primary
      result
    })
    oa_publish_tables(descriptive$tables, descriptive_dir)
    figure_dir <- file.path(output_root, "figures")
    dir.create(figure_dir, showWarnings = FALSE)
    rp_publication_figure(tables$PHQ_profiles, file.path(figure_dir, "paired_phq_profiles.png"))
    comparison <- compare_reference()
    utils::write.csv(comparison$tables, file.path(output_root, "reference_comparison.csv"), row.names = FALSE)
    oa_json(comparison[c("figure_bytes_identical", "figure_scope")], file.path(output_root, "figure_comparison.json"))
    assert(identical(portable_preflight(), identity), "Inputs changed during analysis.")
    files <- list.files(output_root, recursive = TRUE, full.names = TRUE)
    files <- files[!dir.exists(files) & !grepl("/(RUNNING_LOCK/|STATUS\\.json$)", files)]
    relative <- substring(files, nchar(output_root) + 2L)
    oa_json(list(completed = TRUE, identity_sha256 = digest::digest(identity, algo = "sha256"),
      outcome_fits = specs$stage_id, imputation = FALSE, SEM = FALSE, automatic_selection = FALSE,
      elapsed_seconds = proc.time()[["elapsed"]] - started, publication_figure_bytes_identical = comparison$figure_bytes_identical,
      outputs = setNames(as.list(vapply(files, pm_sha, "")), relative)), marker_path)
    status("completed")
    message("Completed: raw-source master, three fixed fits, all current aggregate tables, participant descriptives and Figure 1.")
    invisible(TRUE)
  }, error = function(error) {
    status("failed", conditionMessage(error))
    stop(conditionMessage(error), call. = FALSE)
  }, interrupt = function(error) {
    status("interrupted", "Completed checkpoints preserved. No automatic retry.")
    stop("Interrupted; completed checkpoints preserved.", call. = FALSE)
  })
}
run_portable()
