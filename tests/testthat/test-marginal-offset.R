test_that("plot_marginal_effects() does not draw panels for offset variables", {
  set.seed(3)
  dat <- data.frame(a = rnorm(40), effort = sample(4:12, 40, replace = TRUE))
  dat$y <- stats::rpois(40, exp(0.3 + 0.4 * dat$a) * dat$effort)
  mod <- glm(y ~ a + offset(log(effort)), family = poisson(), data = dat)
  obj <- structure(
    list(opt_mod = mod,
         scl_params = list(mean = c(a = 0), sd = c(a = 1))),
    class = "multiScaleR"
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(obj)
  expect_equal(names(plots), "a")
  expect_equal(
    plots$a$data$fit[1],
    as.numeric(predict(mod,
                       newdata = data.frame(a = min(dat$a),
                                            effort = mean(dat$effort)),
                       type = "response")),
    tolerance = 1e-8
  )
})

test_that("a variable in both a model term and an offset retains its panel", {
  set.seed(4)
  dat <- data.frame(a = rnorm(40), effort = sample(4:12, 40, replace = TRUE))
  dat$y <- stats::rpois(40, exp(0.3 + 0.4 * dat$a + 0.1 * dat$effort) * dat$effort)
  mod <- glm(y ~ a + effort + offset(log(effort)), family = poisson(), data = dat)
  obj <- structure(
    list(opt_mod = mod,
         scl_params = list(mean = c(a = 0), sd = c(a = 1))),
    class = "multiScaleR"
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(obj)
  expect_named(plots, c("a", "effort"))
})
