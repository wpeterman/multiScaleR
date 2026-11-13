source("tools/benchmarks/spatial_field_cv_larger/run.R")
source("tools/benchmarks/spatial_field_cv_fitted/run.R")

# Convert the public fold and declared parameter vector to the reference's
# node-ordered data. There is no dependence on the package's Gaussian mode.
reference_data <- function(d, x, known) {
  ai <- d$train_rows
  ti <- d$test_rows
  nodes <- unique(d$node[c(ai, ti)])
  beta <- if (known) d$beta else c(x$fitted_intercept, x$fitted_slope)
  range <- if (known) d$range else x$fitted_range
  sd <- if (known) d$field_sd else x$fitted_field_sd
  dispersion <- if (known) d$dispersion else x$fitted_dispersion
  eta <- beta[1L] + beta[2L] * d$data$habitat + log(d$data$exposure)
  power <- if (d$covariance == "exponential") 1 else 2
  lower <- sd * t(chol(exp(-(as.matrix(dist(d$coords[nodes, , drop = FALSE]))/range)^power)))
  list(train_y = d$data$y[ai], test_y = d$data$y[ti], train_eta = eta[ai], test_eta = eta[ti],
    train_xy = d$xy[ai, , drop = FALSE], test_xy = d$xy[ti, , drop = FALSE],
    groups = match(d$node[c(ai, ti)], nodes), lower = lower, family = d$family,
    n_train = length(unique(d$node[ai])), range = range, dispersion = dispersion)
}

reference_pit <- function(d, states, subsets, uniforms) {
  train <- seq_len(d$n_train)
  test <- setdiff(seq_len(ncol(d$lower)), train)
  assessment <- length(d$train_y) + seq_along(d$test_y)
  nodes <- d$groups[assessment] - d$n_train
  K <- tcrossprod(d$lower)
  projection <- solve(K[train, train, drop = FALSE], K[train, test, drop = FALSE])
  innovation <- diag(K[test, test, drop = FALSE] -
    t(projection) %*% K[train, test, drop = FALSE])
  moments <- lapply(states[[1L]], function(chain) {
    sweep((chain$field[, train, drop = FALSE] %*% projection)[, nodes, drop = FALSE],
      2L, d$test_eta, "+")
  })
  at_order <- function(order) {
    rule <- normal_rule(order)
    lapply(moments, function(fixed) {
      vapply(seq_along(nodes), function(i) {
        mu <- exp(outer(fixed[, i], sqrt(max(0, innovation[nodes[i]])) * rule$nodes, "+"))
        y <- d$test_y[i]
        if (d$family == "poisson") {
          below <- ppois(y - 1, mu)
          at <- dpois(y, mu)
        } else {
          below <- pnbinom(y - 1, size = d$dispersion, mu = mu)
          at <- dnbinom(y, size = d$dispersion, mu = mu)
        }
        drop((matrix(below, nrow = nrow(mu)) +
          uniforms[i] * matrix(at, nrow = nrow(mu))) %*% rule$weights)
      }, numeric(nrow(fixed)))
    })
  }
  estimate <- function(s, pit) colMeans(do.call(rbind,
    lapply(pit[s$chains], function(x) x[s$rows, , drop = FALSE])))
  current <- at_order(64L)
  point <- estimate(subsets$full, current)
  refinement_gap <- abs(point - estimate(subsets$full, at_order(32L)))
  se <- function(groups)
    apply(vapply(groups, estimate, numeric(length(nodes)), pit = current), 1L, sd)/sqrt(length(groups))
  list(pit = point, mcse = pmax(se(subsets$individual), se(subsets$batches)),
    refinement_gap = refinement_gap)
}

fitted_reference_run <- function(root, out, part, total) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  destination <- file.path(out, "reference-scenarios.csv")
  prediction_file <- file.path(out, "reference-predictions.csv")
  diagnostics_file <- file.path(out, "reference-temperatures.csv")
  if (any(file.exists(c(destination, prediction_file, diagnostics_file))))
    stop("Use a fresh directory to retain reference records.")
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  rows <- do.call(rbind, lapply(file.path(root, c("part1", "part2")), function(path)
    read.csv(file.path(path, "scenarios.csv"), stringsAsFactors = FALSE)))
  predictions <- do.call(rbind, lapply(file.path(root, c("part1", "part2")), function(path)
    read.csv(file.path(path, "predictions.csv"), stringsAsFactors = FALSE)))
  cases <- subset(rows, replicate <= 2L & family %in% c("poisson", "nbinom2"))
  stopifnot(nrow(cases) == 12L)
  cases <- cases[order(cases$id), ]
  cases <- cases[((seq_len(nrow(cases)) - 1L) %% total) + 1L == part, , drop = FALSE]
  temperatures <- seq(0, 1, .25)
  chains <- 3L; draws <- 1600L; warmup <- 350L
  subsets <- monte_carlo_replicates(draws, chains)
  writeLines(c(capture.output(sessionInfo()), "",
    sprintf("Temperatures %s, %s chains, %s warmup, %s draws", paste(temperatures, collapse = ", "),
      chains, warmup, draws), "", capture.output(print(tools::md5sum(c(
      "tools/benchmarks/spatial_field_cv_fitted/reference.R",
      "tools/benchmarks/spatial_field_cv_fitted/run.R", "tools/benchmarks/spatial_field_cv_larger/run.R"))))),
    file.path(out, "reference-session.txt"))
  for (i in seq_len(nrow(cases))) {
    x <- cases[i, ]
    d <- fitted_simulate(x)
    p <- predictions[predictions$id == x$id, ]
    stopifnot(isTRUE(all.equal(p$observed, d$data$y[d$test_rows])))
    set.seed(x$seed + 48913L)
    uniforms <- runif(x$n_assessment)
    for (method in c("generating", "fitted")) {
      known <- method == "generating"
      status <- x[[if (known) "generating_status" else "fit_status"]]
      available <- identical(status, "scored")
      if (!available && !known) {
        result <- list(value = NULL, message = x$fit_message, warnings = x$fit_warnings)
      } else {
        result <- capture(function() {
          r <- reference_data(d, x, known)
          states <- sample_temperatures(r, x$seed + 609034L + if (known) 0L else 11033L,
            temperatures, warmup, draws, chains)
          diag <- do.call(rbind, lapply(states, attr, "diagnostics"))
          append_frame(cbind(id = x$id, method = method, diag), diagnostics_file)
          list(diag = diag, score = reference_score_mc(states, temperatures, subsets),
            moments = reference_moments_mc(r, states, subsets),
            pit = reference_pit(r, states, subsets, uniforms))
        })
      }
      z <- result$value
      overlap <- if (is.null(z)) NA_real_ else min(z$score$overlap)
      rhat <- if (is.null(z)) NA_real_ else max(z$diag$max_rhat)
      mcse <- if (is.null(z)) NA_real_ else z$score$mcse
      score_resolved <- !is.null(z) && rhat <= 1.05 && overlap >= .02 && mcse/x$n_assessment <= .015
      append_frame(data.frame(id = x$id, family = x$family, design = x$design, replicate = x$replicate,
        method = method, public_status = status,
        reference_status = if (!available && !known) "no_fit" else if (is.null(z)) "failed" else
          if (score_resolved) "resolved" else "unresolved",
        n_assessment = x$n_assessment, public_loss = x[[paste0(method, "_loss")]],
        reference_loss = if (is.null(z)) NA_real_ else z$score$loss, reference_loss_mcse = mcse,
        max_rhat = rhat, overlap = overlap, message = result$message, warnings = result$warnings),
        destination)
      if (!is.null(z)) {
        m <- z$moments
        analysis_ok <- z$diag$max_rhat[1L] <= 1.05
        append_frame(data.frame(id = x$id, family = x$family, design = x$design, method = method,
          row = p$row, site = p$site, observed = p$observed, uniform = uniforms,
          reference_pit = z$pit$pit, pit_mcse = z$pit$mcse,
          pit_refinement_gap = z$pit$refinement_gap,
          pit_resolved = analysis_ok & z$pit$mcse <= .015 & z$pit$refinement_gap <= .002,
          reference_mean = m$point[, "mean"], mean_mcse = m$mcse[, "mean"],
          mean_resolved = analysis_ok & m$mcse[, "mean"] <= .03 * abs(m$point[, "mean"]),
          reference_variance = m$point[, "variance"], variance_mcse = m$mcse[, "variance"],
          variance_resolved = analysis_ok & m$mcse[, "variance"] <= .1 * m$point[, "variance"],
          public_mean = p[[paste0(method, "_mean")]],
          public_variance = p[[paste0(method, "_variance")]]), prediction_file)
      }
      cat(x$id, method, "reference completed", !is.null(z), "score resolved", score_resolved, "\n")
    }
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) stop("Supply the completed main results and a fresh output directory.")
  fitted_reference_run(args[1L], args[2L], if (length(args) > 2L) as.integer(args[3L]) else 1L,
    if (length(args) > 3L) as.integer(args[4L]) else 1L)
}
