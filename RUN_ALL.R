# From the repository root: Rscript RUN_ALL.R --run
# In RStudio with this repository as the working directory: source("RUN_ALL.R")
local({
  source_files <- Filter(function(x) is.character(x) && length(x) == 1L,
    lapply(sys.frames(), function(x) x$ofile))
  if (length(source_files)) {
    this_file <- tail(source_files, 1L)[[1L]]
    mode <- "--run"
  } else {
    entry <- grep("^--file=", commandArgs(), value = TRUE)
    if (length(entry) != 1L) stop("Use source() or Rscript to run RUN_ALL.R.", call. = FALSE)
    this_file <- sub("^--file=", "", entry)
    args <- commandArgs(trailingOnly = TRUE)
    mode <- if (length(args)) args else "--run"
  }
  if (length(mode) != 1L || !mode %in% c("--preflight", "--run", "--verify"))
    stop("Choose --preflight, --run or --verify.", call. = FALSE)
  root <- dirname(normalizePath(this_file, winslash = "/", mustWork = TRUE))
  project_library <- file.path(root, ".library")
  if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))
  old_library <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
  on.exit(if (is.na(old_library)) Sys.unsetenv("R_LIBS_USER") else Sys.setenv(R_LIBS_USER = old_library), add = TRUE)
  source(file.path(root, "scripts/public_validation.R"), local = TRUE)
  public_dependencies(root)
  public_check_manifest(root)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (!file.exists(rscript)) stop("Rscript is missing from this R installation.", call. = FALSE)
  for (stage in c("scripts/run.R", "scripts/run_secondary.R", "strengthening/run.R")) {
    message("Executing ", stage, " ", mode)
    status <- system2(rscript, c("--vanilla", shQuote(file.path(root, stage)), mode))
    if (status != 0L) stop(stage, " stopped with exit status ", status, ". Completed checkpoints remain local.", call. = FALSE)
  }
  if (mode != "--preflight") public_compare(root)
  message(if (mode == "--preflight") "Preflight complete. No observed outcome models were fitted." else
    "All three stages completed and matched the fixed manuscript reference values.")
})
