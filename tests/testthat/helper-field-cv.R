# Existing declared Gaussian outcomes, with their physical predictor units
# restored. These are software conformance cases, not new recovery outcomes.
field_cv_fixture <- local({
  common <- NULL; cache <- list(); failure <- list()
  function(covariance = "exponential", kernel = FALSE) {
    if (is.null(common)) {
      original <- native_response_field_fixture("gaussian")
      d <- original$data
      for (nm in c("mean", "rough"))
        d[[nm]] <- d[[nm]] * original$object$scl_params$sd[nm] + original$object$scl_params$mean[nm]
      xy <- terra::crds(original$points); rownames(xy) <- rownames(d)
      planned <- spatial_cv_folds(xy, block_size = 40, k = 2L, crs = terra::crs(original$points), ids = rownames(d))
      folds <- lapply(planned, identity)
      # Exclude one complete training location in the first split. It supplies
      # an explicit buffer-response leakage probe, without using those outcomes.
      position <- as.character(glmmTMB::numFactor(xy))
      first <- match(folds[[1L]]$analysis[1L], rownames(d))
      excluded <- rownames(d)[position == position[first]]
      folds[[1L]]$analysis <- setdiff(folds[[1L]]$analysis, excluded)
      common <<- list(data = d, xy = xy, crs = terra::crs(original$points), folds = folds,
        original = original, excluded = excluded)
    }
    key <- paste(covariance, kernel)
    if (!is.null(failure[[key]])) stop(failure[[key]])
    if (is.null(cache[[key]])) {
      audit <- gdr_test_audit_path(paste0("field-cv-", covariance, "-", kernel, "-"))
      result <- tryCatch(cross_validate_field(y ~ mean + rough + z + offset(log(exposure)),
        common$data, common$xy, common$folds, audit_dir = audit, crs = common$crs,
        covariance = covariance, keep_fits = TRUE,
        kernel_inputs = if (kernel) common$original$inputs else NULL,
        optim_args = if (kernel) list(start_strategy = "single") else list(), max_refits = 200L),
        error = function(e) {failure[[key]] <<- conditionMessage(e); stop(e)})
      if (!result$complete) {
        failure[[key]] <<- paste("Declared field CV fixture failed:",
          paste(vapply(result$errors, function(x) if (is.null(x)) "" else x$message, ""), collapse = "; "))
        stop(failure[[key]])
      }
      cache[[key]] <<- c(common, list(result = result, audit = audit))
    }
    cache[[key]]
  }
})

field_cv_dense_reference <- function(f, name) {
  x <- f$result; s <- x$folds[[name]]; d <- f$data
  ai <- match(s$analysis, rownames(d)); ti <- match(s$assessment, rownames(d))
  model <- x$fits[[name]]; parameters <- x$training[[name]]$parameters
  formula <- y ~ mean + rough + z + offset(log(exposure))
  X <- model.matrix(formula, d); beta <- glmmTMB::fixef(model)$cond
  fixed <- drop(X %*% beta) + log(d$exposure)
  xy <- unname(f$xy)
  distances <- as.matrix(dist(xy))
  power <- if (x$settings$covariance == "exponential") 1 else 2
  K <- parameters$field_SD^2 * exp(-(distances/parameters$range)^power)
  C <- K[ai, ai]; diag(C) <- diag(C) + parameters$observation_SD^2
  cross <- K[ti, ai]
  prediction <- fixed[ti] + drop(cross %*% solve(C, d$y[ai] - fixed[ai]))
  V <- K[ti, ti] - cross %*% solve(C, t(cross))
  diag(V) <- diag(V) + parameters$observation_SD^2
  residual <- d$y[ti] - prediction
  score <- .5 * (length(ti) * log(2 * pi) + as.numeric(determinant(V, logarithm = TRUE)$modulus) +
    drop(crossprod(residual, solve(V, residual))))
  list(prediction = prediction, covariance = V, score = score, rows = ti, training = ai)
}
