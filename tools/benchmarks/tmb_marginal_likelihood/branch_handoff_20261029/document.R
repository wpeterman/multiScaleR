args <- commandArgs(TRUE)
snapshot <- args[[1L]]
study <- args[[2L]]
.libPaths(c(args[[3L]], args[[4L]], .libPaths()))
setwd(snapshot)
paths <- c("NAMESPACE", list.files("man", full.names = TRUE))
before <- tools::md5sum(paths)
roxygen2::roxygenise(".", load_code = function(path) {
  pkgload::load_all(path, helpers = FALSE, compile = TRUE)
  asNamespace("multiScaleR")
})
after <- tools::md5sum(paths)
stopifnot(identical(before, after))
record <- list(passed = TRUE, generated_documentation_unchanged = TRUE,
               libraries = .libPaths(), R_version = R.version.string)
temporary <- file.path(study, "documentation01.json.tmp")
writeLines(jsonlite::toJSON(record, auto_unbox = TRUE), temporary)
stopifnot(file.rename(temporary, file.path(study, "documentation01.json")))
