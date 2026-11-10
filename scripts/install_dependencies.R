# Install pinned packages into this checkout's ignored .library directory.
entry <- grep("^--file=", commandArgs(), value = TRUE)
if (length(entry) != 1L) stop("Run with Rscript scripts/install_dependencies.R.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", entry), winslash = "/")))
if (as.character(getRversion()) != "4.5.1") stop("Install and run R 4.5.1 first.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) && args != "--isolated")) stop("Use no argument or --isolated.")
lib <- file.path(root, ".library"); dir.create(lib, showWarnings = FALSE)
if (identical(args, "--isolated")) .libPaths(c(lib, .Library)) else .libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
m <- read.csv(file.path(root, "scripts/dependency_versions.csv"), stringsAsFactors = FALSE)
order <- c("lattice", "Matrix", "survival", "Rcpp", "RcppArmadillo", "minqa", "numDeriv", "DBI", "mitools", "survey", "jsonlite", "digest", "foreign")
metadata <- list()
available <- function(type) {
  if (is.null(metadata[[type]])) metadata[[type]] <<- available.packages(type = type)
  metadata[[type]]
}
exact_current <- function(name, version, type) {
  a <- available(type)
  name %in% rownames(a) && utils::compareVersion(a[name, "Version"], version) == 0
}
for (name in order) {
  version <- m$version[match(name, m$package)]
  installed <- tryCatch(as.character(packageVersion(name)), error = function(e) "")
  if (nzchar(installed) && utils::compareVersion(installed, version) == 0) {
    message("Already available: ", name, " ", version); next
  }
  # Prefer an exact current binary; never substitute a newer version.
  binary <- if (.Platform$OS.type == "windows") "win.binary" else if (Sys.info()[["sysname"]] == "Darwin") .Platform$pkgType else NULL
  if (!is.null(binary) && exact_current(name, version, binary)) {
    message("Installing exact binary: ", name, " ", version)
    install.packages(name, lib = lib, type = binary, dependencies = FALSE)
  } else {
    current <- exact_current(name, version, "source")
    url <- if (current) paste0("https://cran.r-project.org/src/contrib/", name, "_", version, ".tar.gz") else
      paste0("https://cran.r-project.org/src/contrib/Archive/", name, "/", name, "_", version, ".tar.gz")
    message("Installing pinned source: ", name, " ", version)
    tryCatch(install.packages(url, repos = NULL, type = "source", lib = lib),
      error = function(e) stop(conditionMessage(e), ". Compiled sources require Rtools45 on Windows or a suitable platform compiler."))
  }
  installed <- tryCatch(as.character(packageVersion(name)), error = function(e) "")
  if (!nzchar(installed) || utils::compareVersion(installed, version) != 0)
    stop("Exact installation failed for ", name, ". Compiled sources require Rtools45 on Windows or a platform compiler; also check CRAN archive access. No alternate version was accepted.")
}
source(file.path(root, "scripts/public_validation.R"))
public_dependencies(root)
message("All pinned dependencies verified.")
