# Exact checkpoint helpers only. No observed-analysis driver is loaded.
oa_assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
oa_sha <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
oa_json <- function(x, path) jsonlite::write_json(x, path, pretty = TRUE,
  auto_unbox = TRUE, na = "null", null = "null", digits = NA)

oa_save_once <- function(value, path) {
  oa_assert(!file.exists(path), paste("Refusing to replace saved artifact:", path))
  temporary <- tempfile(pattern = "pending_", tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(value, temporary, version = 3L)
  oa_assert(file.rename(temporary, path), paste("Could not save:", path))
  invisible(path)
}

oa_bind_tables <- function(tables) {
  tables <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, tables)
  if (!length(tables)) return(data.frame())
  columns <- unique(unlist(lapply(tables, names)))
  tables <- lapply(tables, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[columns]
  })
  do.call(rbind, tables)
}

oa_write_csv_once <- function(value, path) {
  oa_assert(!file.exists(path), paste("Refusing to replace saved table:", path))
  utils::write.csv(value, path, row.names = FALSE, na = "")
  invisible(path)
}

oa_stage_checkpoint <- function(path, binding, compute) {
  if (file.exists(path)) {
    saved <- readRDS(path)
    oa_assert(identical(saved$binding, binding) && is.list(saved$result) &&
      is.data.frame(saved$result$table) && nrow(saved$result$table) > 0L,
      "The saved model checkpoint is incompatible. Preserve it and review the changed inputs.")
    return(list(result = saved$result, reused = TRUE))
  }
  result <- compute()
  oa_assert(is.list(result) && is.data.frame(result$table) && nrow(result$table) > 0L,
    "Model returned no contrast table.")
  oa_save_once(list(binding = binding, result = result), path)
  list(result = result, reused = FALSE)
}

oa_acquire_lock <- function(out) {
  parent <- normalizePath(out, winslash = "/", mustWork = TRUE)
  path <- file.path(parent, "RUNNING_LOCK")
  oa_assert(dir.create(path, showWarnings = FALSE),
    "A run lock exists. Do not start a second copy. Check whether the earlier R process is still active.")
  token <- paste(Sys.getpid(), basename(tempfile(pattern = "owner_")), sep = ":")
  owner <- file.path(path, "owner.txt")
  writeLines(token, owner)
  list(parent = parent, path = normalizePath(path, winslash = "/", mustWork = TRUE), token = token)
}

oa_release_lock <- function(lease) {
  # Remove only the exact transient directory created by this invocation.
  # Native R requires recursive=TRUE even when deleting an empty directory.
  if (!dir.exists(lease$path)) return(invisible(TRUE))
  resolved <- normalizePath(lease$path, winslash = "/", mustWork = TRUE)
  parent <- normalizePath(lease$parent, winslash = "/", mustWork = TRUE)
  owner <- file.path(resolved, "owner.txt")
  entries <- list.files(resolved, all.files = TRUE, no.. = TRUE)
  owned <- identical(dirname(resolved), parent) && identical(basename(resolved), "RUNNING_LOCK") &&
    identical(resolved, file.path(parent, "RUNNING_LOCK")) && identical(entries, "owner.txt") &&
    file.exists(owner) && identical(readLines(owner, warn = FALSE), lease$token)
  if (!owned) {
    warning("Run lock ownership or contents changed; it was preserved for inspection.", call. = FALSE)
    return(invisible(FALSE))
  }
  unlink(resolved, recursive = TRUE)
  removed <- !dir.exists(resolved)
  if (!removed) warning("Owned run lock could not be removed; inspect it before restarting.", call. = FALSE)
  invisible(removed)
}
oa_publish_tables <- function(tables, out) {
  # Each table has an immutable RDS receipt, so interrupted publication can resume.
  paths <- character()
  for (name in names(tables)) {
    oa_assert(grepl("^[A-Za-z0-9_]+$", name), "Unsafe output table name.")
    path <- file.path(out, paste0(name, ".csv"))
    receipt <- file.path(out, paste0(name, "_table_receipt.rds"))
    value <- tables[[name]]
    if (file.exists(receipt)) {
      saved <- readRDS(receipt)
      oa_assert(identical(saved$value, value), paste("Saved table content changed:", name))
      if (file.exists(path)) {
        oa_assert(identical(oa_sha(path), saved$sha256), paste("Saved CSV changed:", name))
      } else {
        oa_write_csv_once(value, path)
        oa_assert(identical(oa_sha(path), saved$sha256), paste("Recreated CSV differs:", name))
      }
    } else {
      oa_assert(!file.exists(path), paste("Unreceipted CSV exists; preserve it and review:", path))
      temporary <- tempfile(pattern = "pending_table_", tmpdir = out)
      utils::write.csv(value, temporary, row.names = FALSE, na = "")
      oa_save_once(list(value = value, sha256 = oa_sha(temporary)), receipt)
      oa_assert(file.rename(temporary, path), paste("Could not publish table:", name))
    }
    paths <- c(paths, path, receipt)
  }
  paths
}
