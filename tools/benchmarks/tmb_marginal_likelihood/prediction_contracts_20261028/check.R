Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", MAKEFLAGS = "-j1")
.libPaths(c(normalizePath("tools/benchmarks/joint_scale_assessment/phase31/library"),
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
root <- getwd()
study <- file.path(root, "tools/benchmarks/tmb_marginal_likelihood/prediction_contracts_20261028")
out <- file.path(study, "check01")
stopifnot(!dir.exists(out))
dir.create(out)
validation <- tempfile("msr-prediction-contract-", tmpdir = Sys.getenv("TEMP"))
pkg <- file.path(validation, "multiScaleR")
dir.create(pkg, recursive = TRUE)
files <- unique(c(system2("git", "ls-files", stdout = TRUE),
  system2("git", c("ls-files", "--others", "--exclude-standard", "vignettes", "inst", "src"), stdout = TRUE),
  list.files("R", pattern = "[.]R$", full.names = TRUE),
  list.files("man", pattern = "[.]Rd$", full.names = TRUE),
  list.files("tests", recursive = TRUE, full.names = TRUE)))
files <- files[file.exists(files) & !dir.exists(files) &
  !grepl("^(tools/|[.]git/)|[.](o|dll|so|a)$", files)]
for (file in files) {
  target <- file.path(pkg, file)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(root, file), target))
}
source_hash <- tools::md5sum(files)
stopifnot(identical(unname(source_hash), unname(tools::md5sum(file.path(pkg, files)))))
write.csv(data.frame(file = files, md5 = unname(source_hash)),
  file.path(out, "checked-source-md5.csv"), row.names = FALSE)
writeLines(validation, file.path(out, "check-location.txt"))
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
result <- devtools::check(pkg = pkg, document = FALSE, args = "--no-manual",
  check_dir = file.path(validation, "check"), error_on = "never", quiet = FALSE)
capture.output(result, file = file.path(out, "check-summary.txt"))
logs <- list.files(validation, pattern = "00check.log$", recursive = TRUE, full.names = TRUE)
if (length(logs)) stopifnot(file.copy(logs[1L], file.path(out, "R-CMD-check.log")))
test_logs <- list.files(validation, pattern = "^testthat[.]Rout([.]fail)?$", recursive = TRUE, full.names = TRUE)
if (length(test_logs)) stopifnot(file.copy(test_logs[1L], file.path(out, "installed-tests.log")))
stable <- identical(source_hash, tools::md5sum(files))
record <- list(passed = !length(result$errors) && !length(result$warnings) &&
  !length(result$notes) && stable, errors = result$errors, warnings = result$warnings,
  notes = result$notes, source_stable = stable, checked_source_count = length(files),
  candidate = pkg, check = file.path(validation, "check"), libraries = .libPaths(),
  source_bytes = sum(file.info(file.path(pkg, files))$size), R = R.version.string)
writeLines(jsonlite::toJSON(record, auto_unbox = TRUE, null = "null"), file.path(study, "check01.json.tmp"))
stopifnot(file.rename(file.path(study, "check01.json.tmp"), file.path(study, "check01.json")))
stopifnot(record$passed)
