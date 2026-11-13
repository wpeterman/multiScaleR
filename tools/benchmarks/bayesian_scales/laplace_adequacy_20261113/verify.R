# Verify final archive completeness, fixed sources and independent arithmetic.
root <- "tools/benchmarks/bayesian_scales/laplace_adequacy_20261113"
archive <- file.path(root, "records")
configs <- lapply(file.path(archive, paste0("configuration-", 1:2, ".rds")), readRDS)
stopifnot(identical(configs[[1]], configs[[2]]))
conf <- configs[[1]]
actual <- vapply(names(conf$hash), function(p) digest::digest(file = p, algo = "sha256"), character(1))
changed <- names(actual)[actual != conf$hash]
diagnostic <- "R/bayes_laplace_diagnostic.R"
if (diagnostic %in% changed) {
  snapshot <- file.path(root, "frozen-diagnostic.R")
  old <- digest::digest(file = snapshot, algo = "sha256")
  stopifnot(identical(old, unname(conf$hash[[diagnostic]])))
  code <- function(path) {
    lines <- readLines(path, warn = FALSE)
    lines[!startsWith(lines, "#'")]
  }
  stopifnot(identical(code(snapshot), code(diagnostic)))
  actual[[diagnostic]] <- old
}
for (path in setdiff(changed, diagnostic)) {
  # The review worktree carries other inference updates. Verify original
  # source bytes separately from the package that contains this study.
  snapshot <- file.path(root, "frozen-source", basename(path))
  stopifnot(file.exists(snapshot))
  actual[[path]] <- digest::digest(file = snapshot, algo = "sha256")
}
stopifnot(identical(actual, conf$hash))
ids <- with(conf$plan, expand.grid(case = case, method = conf$methods))
shared <- conf$plan$case[conf$plan$years == 1L]
expected <- with(ids[!(ids$case %in% shared & ids$method == "marginal"), ],
  paste(case, method, "rds", sep = "."))
files <- list.files(archive, pattern = "[.]rds$")
files <- setdiff(files, paste0("configuration-", 1:2, ".rds"))
stopifnot(length(expected) == 95L, setequal(files, expected))
records <- lapply(file.path(archive, files), readRDS)
stopifnot(all(vapply(records, function(x) identical(x$hash, digest::digest(actual)), logical(1))))
stopifnot(identical(files, vapply(records,
  function(x) paste(x$case, x$method, "rds", sep = "."), character(1))))
names(records) <- files
p <- read.csv(file.path(root, "comparisons.csv"))
f <- read.csv(file.path(root, "fits.csv"))
stopifnot(nrow(f) == 105L, !anyDuplicated(f[c("case", "method")]),
  sum(f$status == "outside_scope") == 10L, !any(f$status == "unrun"),
  nrow(p) == 630L, !anyDuplicated(p[c("case", "method", "covariate", "row", "threshold")]),
  all(abs(p$abs_error - abs(p$approximation - p$nuts_estimate)) < 1e-13),
  all(p$fit_status[p$audit_status != "unresolved"] == "eligible"))
key <- c("covariate", "group", "row", "threshold")
match_events <- function(a, b) {
  # The benchmark labels and fixed numeric thresholds cannot contain this
  # separator. Preserve the original ordering when comparing archived rows.
  id <- function(x) do.call(paste, c(unname(as.list(x[key])), list(sep = "\\r")))
  stopifnot(!anyDuplicated(id(a)), !anyDuplicated(id(b)), setequal(id(a), id(b)))
  match(id(a), id(b))
}
for (record in records) {
  fit <- f[f$case == record$case & f$method == record$method, , drop = FALSE]
  stopifnot(nrow(fit) == 1L, identical(fit$status, record$status))
  if (record$method == "nuts") next
  reference <- records[[paste(record$case, "nuts", "rds", sep = ".")]]
  saved <- record$audit$probabilities
  original <- reference$probability
  actual <- p[p$case == record$case & p$method == record$method, , drop = FALSE]
  stopifnot(nrow(actual) == nrow(saved), nrow(actual) == nrow(original))
  names(saved)[names(saved) == "status"] <- "audit_status"
  names(original)[match(c("estimate", "mcse"), names(original))] <-
    c("nuts_estimate", "nuts_mcse")
  for (table in list(saved, original[c(key, "nuts_estimate", "nuts_mcse")])) {
    fields <- names(table)
    archived <- table[match_events(actual, table), , drop = FALSE]
    stopifnot(isTRUE(all.equal(actual[fields], archived, tolerance = 1e-13,
      check.attributes = FALSE)))
  }
  stopifnot(all(actual$fit_status == record$status),
    all(actual$nuts_status == reference$status))
  n <- 8000L
  z <- 3
  a <- actual$approximation
  center <- (a + z^2 / (2 * n)) / (1 + z^2 / n)
  radius <- z * sqrt(a * (1 - a) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  lo <- center - radius - actual$nuts_estimate - z * actual$nuts_mcse
  hi <- center + radius - actual$nuts_estimate + z * actual$nuts_mcse
  oracle <- rep("unresolved", nrow(actual))
  eligible <- record$status == "eligible" && reference$status == "eligible"
  if (eligible) {
    oracle[which(lo > .05 | hi < -.05)] <- "discrepancy"
    oracle[which(lo >= -.05 & hi <= .05)] <- "within_tolerance"
  }
  stopifnot(identical(actual$oracle_status, oracle))
}
stopifnot(all(p$reference_abs_error[is.finite(p$reference_abs_error)] >= 0),
  all(p$indicator_ess[p$audit_status != "unresolved"] >= 400),
  all(p$split_reference_difference[p$audit_status != "unresolved"] <= .05))
cat(sprintf("Verified %d planned fits and %d paired event rows.\n", length(expected), nrow(p)))
