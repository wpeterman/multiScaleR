# Bounded native-constructor fixtures, not a recovery or simulation study. Fixed
# fits reuse the declared C fixture's data. The mixed fixture uses one fixed seed
# and explicit balanced group effects; failed fits are retained, never reseeded.
gdr_tmb_refit_fixture <- local({
  cache <- list()
  function(kind = "P") {
    if (!is.null(cache[[kind]])) return(cache[[kind]])
    testthat::skip_if(as.character(utils::packageVersion("unmarked")) != "1.5.2" ||
      as.character(utils::packageVersion("TMB")) != "1.9.25", "Native TMB provider needs qualification")
    if (kind != "mixed") {
      original <- gdr_refit_fixture(kind)$model
      frame <- original@data
      starts <- unname(unmarked::coef(original))
      if (kind == "NB") starts <- starts[c(1, 2, 4, 5, 3)]
      forms <- original@formlist
      mixture <- kind
    } else {
      existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      saved <- if (existed) get(".Random.seed", .GlobalEnv) else NULL
      on.exit(if (existed) assign(".Random.seed", saved, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
      set.seed(6127L)
      n <- 160L
      group <- factor(rep(seq_len(8L), each = 20L))
      observer <- factor(rep(seq_len(4L), length.out = n))
      hab <- rep(seq(-1.5, 1.5, length.out = 20L), 8L)
      state <- rep(c(-.9, -.65, -.4, -.1, .1, .4, .65, .9), each = 20L)
      detection <- rep(c(-.3, -.1, .1, .3), length.out = n)
      abundance <- rpois(n, 900 * exp(.3 * hab + state) * pi * .075^2)
      breaks <- c(0, 25, 50, 75); periods <- c(1, 2, 1)
      removal <- .7^c(0, 1, 3) * (1 - .7^periods)
      yd <- yr <- matrix(0, n, 3L)
      for (i in seq_len(n)) {
        sigma <- 40 * exp(detection[i])
        distance <- 2 * sigma^2 / 75^2 * diff(-exp(-breaks^2 / (2 * sigma^2)))
        probability <- as.vector(outer(distance, removal))
        joint <- matrix(rmultinom(1L, abundance[i], c(probability, 1 - sum(probability)))[1:9], 3L)
        yd[i, ] <- rowSums(joint); yr[i, ] <- colSums(joint)
      }
      frame <- unmarked::unmarkedFrameGDR(yd, yr, dist.breaks = breaks,
        unitsIn = "m", period.lengths = periods,
        siteCovs = data.frame(hab = hab, block = group, observer = observer))
      forms <- lapply(c(lambda = "~hab + (1 | block)", phi = "~1",
        dist = "~1 + (1 | observer)", rem = "~1"), stats::as.formula, env = baseenv())
      starts <- c(log(900), .3, log(40), qlogis(.3), log(.6), log(.25))
      mixture <- "P"
    }
    settings <- list(method = "BFGS", se = TRUE, threads = 1L,
      control = list(maxit = 2500L, reltol = 1e-10))
    args <- list(lambdaformula = forms$lambda, phiformula = forms$phi,
      distanceformula = forms$dist, removalformula = forms$rem, data = frame,
      K = 400, keyfun = "halfnorm", mixture = mixture, output = "density",
      unitsOut = "kmsq", engine = "TMB", threads = 1L, starts = starts, method = settings$method)
    record <- paste0(gdr_test_audit_path(paste0("initial-TMB-", kind, "-")), ".rds")
    model <- tryCatch(.msr_gdr_TMB_refit(args, list(dots = list(control = settings$control))),
      error = function(e) {
        .msr_gdr_atomic_audit(list(kind = kind, status = "initial_constructor_failed",
          error = conditionMessage(e)), record)
        stop(e)
      })
    .msr_gdr_atomic_audit(list(kind = kind, convergence = model@opt$convergence,
      logLik = as.numeric(unmarked::logLik(model)),
      diagnostic = attr(model, "msr_gdr_tmb_diagnostics")), record)
    snapshot <- gdr_start_values(model)
    cold <- snapshot$starts
    names(cold) <- paste0(snapshot$map$block, "[", snapshot$map$within_block, "]")
    cold[snapshot$map$block == "beta_lambda"] <- c(log(700), 0)
    cold[snapshot$map$block == "beta_dist"] <- log(30)
    cold[snapshot$map$block == "beta_rem"] <- qlogis(.5)
    cold[startsWith(snapshot$map$block, "lsigma_")] <- log(.4)
    if (kind == "NB") cold[snapshot$map$block == "beta_alpha"] <- log(1.5)
    updated <- frame; updated@siteCovs$hab <- .8 * frame@siteCovs$hab + .1
    cache[[kind]] <<- list(model = model, data = updated, settings = settings,
      starts = list(cold = cold), constructor_args = args)
    cache[[kind]]
  }
})
