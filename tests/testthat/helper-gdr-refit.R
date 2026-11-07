# Three declared single-primary fixtures, one fixed seed per generating mixture.
# This test generator is not an extension of the public Poisson-only simulator.
gdr_test_audit_path <- function(prefix) {
  root <- Sys.getenv("MSR_GDR_TEST_AUDIT_ROOT", unset = "")
  if (!nzchar(root)) return(tempfile(prefix))
  if (!dir.exists(root)) dir.create(root, recursive = TRUE)
  tempfile(prefix, tmpdir = root)
}

gdr_contract_model <- local({
  model <- NULL
  function() {
    if (is.null(model)) model <<- unmarked::gdistremoval(lambdaformula = ~hab,
      data = gdr_guard_fixture()$P@data, output = "abund", K = 40,
      starts = c(log(8), 0, log(40), qlogis(0.3)), engine = "C",
      threads = 1, se = FALSE, control = list(maxit = 1))
    model
  }
})

gdr_refit_fixture <- local({
  cache <- list()
  function(mixture) {
    if (!is.null(cache[[mixture]])) return(cache[[mixture]])
    existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
    old_seed <- if (existed) get(".Random.seed", .GlobalEnv) else NULL
    on.exit(if (existed) assign(".Random.seed", old_seed, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
    set.seed(switch(mixture, P = 6041L, NB = 6042L, ZIP = 6043L))
    n <- 80L
    x <- seq(-1.5, 1.5, length.out = n)
    area <- pi * 75^2 / 1e6
    mean_n <- 800 * exp(0.3 * x) * area
    abundance <- if (mixture == "NB") rnbinom(n, mu = mean_n, size = 2) else rpois(n, mean_n)
    if (mixture == "ZIP") abundance <- abundance * rbinom(n, 1, 0.6)
    breaks <- c(0, 25, 50, 75)
    durations <- c(1, 2, 1)
    distance <- 2 * 40^2 / 75^2 * diff(-exp(-breaks^2 / (2 * 40^2)))
    before <- c(0, utils::head(cumsum(durations), -1))
    removal <- 0.7^before * (1 - 0.7^durations)
    observed <- as.vector(outer(distance, removal))
    probabilities <- c(observed, 1 - sum(observed))
    yd <- yr <- matrix(0, n, 3)
    for (i in seq_len(n)) {
      joint <- matrix(rmultinom(1, abundance[i], probabilities)[seq_len(9)], 3, 3)
      yd[i, ] <- rowSums(joint)
      yr[i, ] <- colSums(joint)
    }
    frame <- unmarked::unmarkedFrameGDR(yd, yr, dist.breaks = breaks, unitsIn = "m",
      period.lengths = durations, siteCovs = data.frame(hab = x))
    alpha <- switch(mixture, P = numeric(), NB = log(2), ZIP = qlogis(0.4))
    # Explicit formula environments avoid serializing a test helper's package scope.
    one <- stats::as.formula("~1", env = baseenv())
    lambda <- stats::as.formula("~hab", env = baseenv())
    warnings <- character()
    initial_record <- paste0(gdr_test_audit_path(paste0("initial-", mixture, "-")), ".rds")
    model <- tryCatch(withCallingHandlers(unmarked::gdistremoval(lambdaformula = lambda, phiformula = one,
      distanceformula = one, removalformula = one, data = frame,
      output = "density", unitsOut = "kmsq", K = 400, mixture = mixture,
      starts = c(log(800), 0.3, alpha, log(40), qlogis(0.3)),
      engine = "C", threads = 1, method = "BFGS", se = TRUE,
      control = list(maxit = 2500, reltol = 1e-10)),
      warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
      error = function(e) {
        .msr_gdr_atomic_audit(list(mixture = mixture, status = "initial_constructor_failed", warnings = warnings,
          error = conditionMessage(e), seed = switch(mixture, P = 6041L, NB = 6042L, ZIP = 6043L)), initial_record)
        stop(e)
      })
    .msr_gdr_atomic_audit(list(mixture = mixture, status = "initial_constructor_returned", warnings = warnings,
      convergence = model@opt$convergence, logLik = as.numeric(unmarked::logLik(model)),
      coefficients = unmarked::coef(model), covariance = .msr_gdr_variance_status(model, TRUE),
      seed = switch(mixture, P = 6041L, NB = 6042L, ZIP = 6043L)), initial_record)
    cold <- unmarked::coef(model)
    cold[] <- 0
    cold["lambda(Int)"] <- log(600)
    cold["dist(Int)"] <- log(25)
    cold["rem(Int)"] <- qlogis(0.6)
    if (mixture != "P") cold["alpha(alpha)"] <- if (mixture == "NB") log(1.5) else qlogis(0.2)
    other <- cold
    other["lambda(Int)"] <- log(1000)
    other["lambda(hab)"] <- 0.5
    other["dist(Int)"] <- log(60)
    other["rem(Int)"] <- qlogis(0.15)
    updated <- frame
    updated@siteCovs$hab <- 0.8 * x + 0.1
    cache[[mixture]] <<- list(model = model, data = updated,
      starts = list(cold = cold, other = other), seed = switch(mixture, P = 6041L, NB = 6042L, ZIP = 6043L))
    cache[[mixture]]
  }
})

# One actual objective callback on one worker, with the same source/archive
# namespace as the driver. Return compact summaries, not a fitted-model array.
gdr_refit_worker_probe <- function(callback, model, data) {
  cluster <- parallel::makePSOCKcluster(1L)
  on.exit(parallel::stopCluster(cluster))
  path <- getNamespaceInfo(asNamespace("multiScaleR"), "path")
  source <- !file.exists(file.path(path, "R", "multiScaleR.rdb"))
  setup <- function(libraries, source, path) {
    .libPaths(libraries)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    if (source) pkgload::load_all(path, quiet = TRUE) else
      library("multiScaleR", lib.loc = dirname(path), character.only = TRUE)
    TRUE
  }
  environment(setup) <- baseenv()
  parallel::clusterCall(cluster, setup, .libPaths(), source, path)
  job <- function(callback, model, data) {
    fit <- callback(model, data, list(return_model = FALSE, covs = "hab"))
    list(logLik = as.numeric(unmarked::logLik(fit)),
      K = get(".msr_parameter_count", asNamespace("multiScaleR"))(fit),
      audit_file = attr(fit, "msr_gdr_refit_audit")$path)
  }
  environment(job) <- baseenv()
  parallel::clusterCall(cluster, job, callback, model, data)[[1L]]
}
