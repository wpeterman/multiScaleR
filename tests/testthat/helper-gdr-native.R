# Three fixed native C constructors only; maxit=0 supplies likelihood oracles,
# not converged fits or evidence about recovery. No TMB optimizer is launched.
gdr_guard_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    yd <- rbind(c(0, 0), c(1, 0), c(0, 2), c(2, 1), c(0, 0))
    yr <- rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 1), c(1, 1, 1), c(0, 0, 0))
    frame <- unmarked::unmarkedFrameGDR(yd, yr, dist.breaks = c(0, 25, 75),
      unitsIn = "m", period.lengths = c(1, 2, 1),
      siteCovs = data.frame(hab = seq(-1, 1, length.out = 5), exposure = rep(0.4, 5)))
    cache <<- stats::setNames(lapply(c("P", "NB", "ZIP"), function(mixture) {
      alpha <- switch(mixture, P = numeric(), NB = log(2), ZIP = qlogis(0.2))
      unmarked::gdistremoval(lambdaformula = ~hab, data = frame, K = 40,
        output = "abund", mixture = mixture, engine = "C", se = FALSE,
        starts = c(log(8), 0, alpha, log(40), qlogis(0.3)), control = list(maxit = 0))
    }), c("P", "NB", "ZIP"))
    cache
  }
})


# Serializing formula/nll environments traverses the helper's package scope
# and emits loader warnings. Snapshot fitted slots without serializing them;
# direct identity still checks data, coefficients, covariance, optimizer and
# settings, plus the original formula and likelihood-function objects.
gdr_guard_state <- function(model) {
  slots <- methods::slotNames(model)
  stats::setNames(lapply(slots, function(name) methods::slot(model, name)), slots)
}
