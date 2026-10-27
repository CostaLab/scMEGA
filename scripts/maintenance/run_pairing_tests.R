# Focused source tests are explicitly separate from package installation/check.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
Sys.setenv(SCMEGA_SOURCE_TESTS = "true", SCMEGA_SOURCE_ROOT = root)
stopifnot(requireNamespace("testthat"), requireNamespace("SeuratObject"))
if (Sys.getenv("SCMEGA_REQUIRE_INTEGRATION") == "true") {
  stopifnot(requireNamespace("Signac"), requireNamespace("GenomicRanges"))
}
testthat::set_max_fails(Inf)
test.root <- normalizePath(if (length(args) > 1) args[[2]] else root, mustWork = TRUE)
results <- testthat::test_dir(file.path(test.root, "tests", "testthat"), reporter = "summary",
  stop_on_failure = FALSE, filter = "paired")
summary <- as.data.frame(results)
print(summary[, c("test", "failed", "skipped", "error", "warning", "passed")], row.names = FALSE)
if (any(summary$failed > 0 | summary$error)) stop("Pairing regression tests failed.")
if (Sys.getenv("SCMEGA_REQUIRE_INTEGRATION") == "true" && any(summary$skipped)) {
  stop("Required integration tests were skipped.")
}
