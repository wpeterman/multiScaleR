# cv_tab() and select_cv() route any list containing a mixed-model CV result
# to the mixed-model comparison. A candidate of another type used to trigger
# "Every fold must be scored", which named the wrong cause.

test_that("mixing mixed-model and GLM CV results names the incompatible candidate", {
  skip_if_not_installed("glmmTMB")
  set.seed(10)
  site <- factor(rep(sprintf("site%02d", 1:12), each = 6))
  surveys <- data.frame(site, habitat = rnorm(length(site)))
  surveys$count <- rpois(nrow(surveys),
    exp(0.9 + 0.4 * surveys$habitat + rnorm(nlevels(site), sd = 0.7)[site]))
  folds <- grouped_cv_folds(surveys$site, k = 3, ids = rownames(surveys))

  mixed <- suppressWarnings(cross_validate_mixed(count ~ habitat + (1 | site),
    surveys, folds, family = "poisson", prediction = "new_group"))
  glm_cv <- cross_validate_glm(count ~ habitat, surveys, folds, family = "poisson")
  expect_true(mixed$complete)
  expect_true(glm_cv$complete)

  # Both complete, so the message must name the type mismatch and the candidate
  expect_error(cv_tab(list(mixed = mixed, glm = glm_cv), metric = "rmse"),
    "only with each other. Not a mixed-model CV result: glm.", fixed = TRUE)
  expect_error(select_cv(list(mixed = mixed, glm = glm_cv), criterion = "AIC"),
    "Not a mixed-model CV result: glm.", fixed = TRUE)

  # An incomplete mixed-model run still reports the fold problem
  failed <- mixed
  failed$complete <- FALSE
  expect_error(cv_tab(list(a = mixed, b = failed), metric = "rmse"),
    "Every fold must be scored")
})
