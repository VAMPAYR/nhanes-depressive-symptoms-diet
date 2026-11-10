# Download only the official public files listed in the root manifest.
entry <- grep("^--file=", commandArgs(), value = TRUE)
if (length(entry) != 1L) stop("Run with Rscript scripts/download_data.R --download or --verify.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", entry), winslash = "/")))
if (dir.exists(file.path(root, ".library"))) .libPaths(c(file.path(root, ".library"), .libPaths()))
if (!requireNamespace("digest", quietly = TRUE)) stop("Run scripts/install_dependencies.R first.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args %in% c("--download", "--verify")) stop("Choose --download or --verify.")
m <- read.csv(file.path(root, "data_manifest.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(m) == 32L, !anyDuplicated(m$path), all(grepl("^https://wwwn[.]cdc[.]gov/", m$url)))
if (any(grepl("(^[/\\\\]|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))", m$path))) stop("Unsafe input manifest path.")
options(timeout = max(600, getOption("timeout", 60)))
for (i in seq_len(nrow(m))) {
  p <- file.path(root, "data", m$path[i])
  if (!file.exists(p)) {
    if (args == "--verify") stop("Missing input: ", m$path[i])
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(pattern = "download_", tmpdir = dirname(p))
    tryCatch({
      status <- download.file(m$url[i], temporary, mode = "wb", method = "libcurl", quiet = FALSE)
      if (status != 0L || !identical(digest::digest(file = temporary, algo = "sha256", serialize = FALSE), m$sha256[i]))
        stop("Download failed or differed from the fixed public input: ", m$path[i], ". HTML/error responses are rejected.")
      if (!file.rename(temporary, p)) stop("Could not place verified input: ", m$path[i])
    }, finally = if (file.exists(temporary)) unlink(temporary))
  }
  if (!identical(digest::digest(file = p, algo = "sha256", serialize = FALSE), m$sha256[i]))
    stop("Existing input differs from manifest; preserve it and inspect: ", m$path[i])
  message("Verified ", i, "/", nrow(m), ": ", m$path[i])
}
message("All 32 official input files verified.")
