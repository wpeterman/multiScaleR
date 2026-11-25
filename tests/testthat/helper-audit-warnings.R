# Cross-validation functions warn when they return unscored folds. Fixtures
# that build a deliberately incomplete result use this to silence that one
# warning; every other warning still reaches the test.
msr_quiet_incomplete <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("returned an incomplete result", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}
