# Four bounded constructions, fixed seeds and retained failures. These fixtures
# qualify prediction recipes, not recovery, new scale fits or a timing claim.
gdr_trained_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function() {
    if (exists("failure", cache, inherits = FALSE)) stop(cache$failure)
    if (exists("value", cache, inherits = FALSE)) return(cache$value)
    cache$value <- tryCatch({
      existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      saved <- if (existed) get(".Random.seed", .GlobalEnv) else NULL
      on.exit(if (existed) assign(".Random.seed", saved, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
      set.seed(6831L)
      n <- 160L; primary <- 2L; durations <- c(1, 2, 1)
      hab <- seq(-1.5, 1.5, length.out = n)
      cover <- factor(rep(c("None", "Forest"), length.out = n))
      contrasts(cover) <- stats::contr.sum(2L)
      wind <- runif(n * primary, -1.5, 1.5)
      temperature <- runif(n * primary, -1, 1)
      effort <- runif(n * primary * 3L, -.8, .8)
      sw <- as.numeric(scale(wind)); se <- as.numeric(scale(effort))
      abundance <- rpois(n, 1000 * exp(.2 * hab + .15 * (cover == "Forest")) * pi * .075^2)
      yd <- matrix(0, n, 4L); yr <- matrix(0, n, 6L)
      for (i in seq_len(n)) for (j in seq_len(primary)) {
        row <- (i - 1L) * primary + j
        sigma <- 40 * exp(.1 * temperature[row] - .05 * temperature[row]^2)
        distance <- 2 * sigma^2 / 75^2 * diff(-exp(-c(0, 25, 75)^2 / (2 * sigma^2)))
        index <- (row - 1L) * 3L + seq_len(3L)
        interval <- 1 - (1 - plogis(-.7 + .1 * se[index]^2))^durations
        removal <- interval * c(1, utils::head(cumprod(1 - interval), -1L))
        probabilities <- as.vector(outer(distance, removal)) * plogis(.3 + .15 * sw[row]^2 + .1 * (wind[row] >= 0))
        joint <- matrix(rmultinom(1L, abundance[i], c(probabilities, 1 - sum(probabilities)))[1:6], 2L)
        yd[i, (j - 1L) * 2L + seq_len(2L)] <- rowSums(joint)
        yr[i, (j - 1L) * 3L + seq_len(3L)] <- colSums(joint)
      }
      frame <- unmarked::unmarkedFrameGDR(yd, yr, numPrimary = primary,
        dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = durations,
        siteCovs = data.frame(hab = hab, `land cover` = cover, check.names = FALSE),
        yearlySiteCovs = data.frame(wind = wind, temperature = temperature, wind_good = wind >= 0),
        obsCovs = data.frame(effort = effort))
      forms <- list(lambda = ~splines::ns(hab, df = 4) * `land cover`,
        phi = ~I(base::scale(wind)^2) + wind_good, dist = ~stats::poly(temperature, 2), rem = ~I(base::scale(effort)^2))
      # The fixture's constructor workspace contains other models and audit
      # helpers. Its formulas need only data and qualified standard functions.
      for (name in names(forms)) environment(forms[[name]]) <- baseenv()
      design <- unmarked:::getDesign(frame, forms)
      starts <- unlist(lapply(c("lambda", "phi", "dist", "rem"), function(type) {
        x <- numeric(ncol(design[[paste0("X_", type)]]))
        x[1L] <- switch(type, lambda = log(1000), phi = qlogis(.65), dist = log(40), rem = qlogis(.3))
        x
      }), use.names = FALSE)
      settings <- list(dots = list(control = list(maxit = 2500L, reltol = 1e-10)))
      args <- list(lambdaformula = forms$lambda, phiformula = forms$phi,
        distanceformula = forms$dist, removalformula = forms$rem, data = frame,
        K = 200L, keyfun = "halfnorm", mixture = "P", output = "density",
        unitsOut = "kmsq", engine = "TMB", threads = 1L, starts = starts, method = "BFGS", se = TRUE)
      models <- list()
      record <- function(name, expression) {
        warnings <- character(); path <- paste0(gdr_test_audit_path(paste0("trained-", name, "-")), ".rds")
        value <- tryCatch(withCallingHandlers(force(expression), warning = function(w)
          warnings <<- c(warnings, conditionMessage(w))), error = function(e) {
            .msr_gdr_atomic_audit(list(name = name, error = conditionMessage(e), warnings = warnings), path)
            stop(e)
          })
        diagnostic <- .msr_gdr_variance_status(value, TRUE)
        .msr_gdr_atomic_audit(list(name = name, warnings = warnings, convergence = value@opt$convergence,
          logLik = as.numeric(unmarked::logLik(value)), covariance = diagnostic,
          parameter_count = .msr_parameter_count(value)), path)
        stopifnot(!length(warnings), as.integer(value@opt$convergence) == 0L,
          identical(diagnostic$status, "available"))
        value
      }
      models$TMB <- record("fixed-TMB", .msr_gdr_TMB_refit(args, settings))
      args$engine <- "C"; args$starts <- gdr_start_values(models$TMB)$starts
      models$C <- record("fixed-C", do.call(unmarked::gdistremoval, c(args, settings$dots)))
      original <- gdr_tmb_refit_fixture("mixed")$model
      mixed_frame <- original@data
      mixed_frame@siteCovs$temperature <- seq(-1, 1, length.out = nrow(mixed_frame@siteCovs))
      mixed_forms <- original@formlist
      mixed_forms$lambda <- ~splines::ns(hab, df = 3) + (1 | block)
      mixed_forms$dist <- ~stats::poly(temperature, 2) + (1 | observer)
      environment(mixed_forms$lambda) <- environment(mixed_forms$dist) <- baseenv()
      old_design <- unmarked:::getDesign(original@data, original@formlist)
      new_design <- unmarked:::getDesign(mixed_frame, mixed_forms)
      old_values <- gdr_start_values(original)
      coefficient <- function(type) unmarked::coef(original, type = type, altNames = FALSE)[original@estimates@estimates[[type]]@fixed]
      beta_lambda <- as.vector(qr.solve(new_design$X_lambda, old_design$X_lambda %*% coefficient("lambda")))
      beta_dist <- c(coefficient("dist")[1L], 0, 0)
      sigma <- old_values$starts[old_values$map$kind == "log_sd"]
      mixed_args <- list(lambdaformula = mixed_forms$lambda, phiformula = mixed_forms$phi,
        distanceformula = mixed_forms$dist, removalformula = mixed_forms$rem, data = mixed_frame,
        K = original@K, keyfun = original@keyfun, mixture = original@mixture,
        output = original@output, unitsOut = original@unitsOut, engine = "TMB", threads = 1L,
        starts = c(beta_lambda, beta_dist, coefficient("rem"), sigma), method = "BFGS", se = TRUE)
      models$mixed <- record("mixed-ns3", .msr_gdr_TMB_refit(mixed_args, settings))
      list(models = models, provider = body(unmarked:::make_mod_matrix),
        snapshots = lapply(models, function(m) list(estimates = m@estimates, opt = m@opt,
          formlist = m@formlist, data = m@data, K = m@K, AIC = m@AIC)))
    }, error = function(e) {cache$failure <- e; stop(e)})
    cache$value
  }
})

# Independent base-R recipes for these particular fixture formulas. Prediction
# covariance/chunking stays with the actual native provider in both routes.
gdr_trained_reference <- function(model, type, rows, backTransform = TRUE, re.form = NULL) {
  raw <- unmarked:::get_covariates(model, type)
  mixed <- !is.null(model@TMB) && length(model@TMB$env$random) > 0L
  if (type == "lambda") {
    basis <- splines::ns(raw$hab, df = if (mixed) 3L else 4L)
    values <- as.data.frame(predict(basis, rows$hab))
    if (mixed) X <- cbind(1, as.matrix(values)) else {
      values$cover <- factor(rows[["land cover"]], levels = levels(raw[["land cover"]]))
      contrasts(values$cover) <- contrasts(raw[["land cover"]])
      X <- model.matrix(~(.) * cover, values)
    }
  } else if (type == "dist") X <- cbind(1, predict(poly(raw$temperature, 2), rows$temperature)) else if (mixed && type == "rem") {
    X <- matrix(1, nrow(rows), 1L)
  } else {
    name <- if (type == "phi") "wind" else "effort"
    X <- cbind(1, ((rows[[name]] - mean(raw[[name]])) / sd(raw[[name]]))^2)
    if (type == "phi") X <- cbind(X, as.integer(rows$wind_good))
  }
  if (is.null(re.form) && unmarked:::has_random(model@formlist[[type]]))
    X <- cbind(X, unmarked:::get_Z(model@formlist[[type]], raw, rows))
  unmarked:::predict_by_chunk(model, type, .95, X, NULL, 70L, backTransform, re.form)
}
