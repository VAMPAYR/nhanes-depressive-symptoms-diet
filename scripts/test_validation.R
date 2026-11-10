# Synthetic fixtures only; real reference values and results are never modified.
entry <- grep("^--file=", commandArgs(), value = TRUE)
if (length(entry) != 1L) stop("Run with Rscript scripts/test_validation.R.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", entry), winslash = "/")))
if (dir.exists(file.path(root, ".library"))) .libPaths(c(file.path(root, ".library"), .libPaths()))
source(file.path(root, "scripts/public_validation.R"))
fixture <- tempfile("nhanes_reference_test_")
dir.create(fixture)
expected <- data.frame(label = c("a", "b"), estimate = c(1, 2))
write.csv(expected, file.path(fixture, "expected.csv"), row.names = FALSE)
write.csv(data.frame(stage = "test", reference = "expected.csv", output = "actual.csv", type = "csv",
  sha256 = public_sha(file.path(fixture, "expected.csv"))), file.path(fixture, "reference_manifest.csv"), row.names = FALSE)
must_fail <- function(label) {
  error <- tryCatch({public_compare(fixture, "test", FALSE); NULL}, error = function(e) e)
  public_assert(inherits(error, "error"), paste("Invalid fixture was accepted:", label))
  message("Expected failure verified: ", label)
}
must_fail("missing required output")
bad <- expected; bad$estimate[2] <- 3
write.csv(bad, file.path(fixture, "actual.csv"), row.names = FALSE)
must_fail("changed numerical result")
bad <- expected; names(bad)[2] <- "renamed_estimate"
write.csv(bad, file.path(fixture, "actual.csv"), row.names = FALSE)
must_fail("changed column schema")
write.csv(expected, file.path(fixture, "actual.csv"), row.names = FALSE)
public_compare(fixture, "test", FALSE)
message("All reference-verifier fixtures passed. Real analysis files were untouched.")
