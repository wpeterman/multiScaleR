raw_binomial_fixture <- function(grouped = FALSE, factor_response = FALSE) {
  n <- 80L
  forest <- rep(seq(.04, .96, length.out = 20L), each = 4L)
  developed <- rep(.05 + .9 * (1 + sin(seq(0, 5, length.out = 20L))) / 2, each = 4L)
  d <- data.frame(y = rep(c(0, 1), n / 2), forest = forest, developed = developed,
    exposure = exp(.1 * forest), cover = factor(rep(c("a", "a", "b", "b"), n / 4)))
  rownames(d) <- sprintf("point-%03d", seq_len(n))
  d$frequency <- if (grouped) rep(c(1, 2), n / 2) else rep(1, n)
  if (grouped) { d$success <- rep(c(1, 2), n / 2); d$failure <- 3 - d$success }
  if (factor_response) d$y <- factor(ifelse(d$y == 1, "yes", "no"), levels = c("no", "yes"))
  formula <- if (grouped) cbind(success, failure) ~
    splines::ns(forest, knots = .5, Boundary.knots = c(0, 1)) +
    splines::ns(developed, knots = .5, Boundary.knots = c(0, 1)) + offset(log(exposure)) else y ~
    splines::ns(forest, knots = .5, Boundary.knots = c(0, 1)) +
    splines::ns(developed, knots = .5, Boundary.knots = c(0, 1)) + offset(log(exposure))
  fit <- glm(formula, data = d, weights = frequency, family = binomial(),
    control = glm.control(epsilon = 1e-11, maxit = 100L))
  reference <- as.matrix(d[c("forest", "developed")])
  context <- .msr_prepare_raw_binomial(fit, colnames(reference), reference, paste(rep("a", 64), collapse = ""),
    c("sigma_forest", "sigma_developed"))
  list(data = d, formula = formula, fit = fit, reference = reference, context = context)
}
