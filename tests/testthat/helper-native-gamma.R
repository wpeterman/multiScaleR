# Deterministic positive-response software fixture; selected scales and the
# original mean/SQ/fixed-SA reference transform are held fixed, not reestimated.
native_gamma_fixture <- local({
  cache <- list()
  function(kind = "glm") {
    if (!is.null(cache[[kind]])) return(cache[[kind]])
    f <- refit_spatial_fixture()
    eta <- 5 + .03 * f$data$mean + .02 * f$data$rough + .01 * f$data$fixed + .02 * f$data$z
    stopifnot(all(is.finite(eta)), all(eta > 0))
    f$data$y <- exp(.01 * sin(seq_along(eta))) / eta
    f$object$opt_mod <- if (kind == "glm")
      glm(y ~ mean + rough + fixed + z, data = f$data, family = Gamma()) else
        spaMM::fitme(y ~ mean + rough + fixed + z, data = f$data,
          family = Gamma(), method = "ML")
    cache[[kind]] <<- f
    f
  }
})
