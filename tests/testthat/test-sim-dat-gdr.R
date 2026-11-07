sim_gdr_fixture <- function(...) {
  args <- list(density = c(100, 200, 0), distance_scale = 40, removal = 0.3,
    dist_breaks = c(0, 25, 75), interval_lengths = c(1, 2, 1),
    interval_units = "minutes", seed = 891L)
  extra <- list(...); args[names(extra)] <- extra
  do.call(sim_dat_gdr, args)
}

test_that("GDR simulation retains paired observations, latent truth, and provenance", {
  sites <- data.frame(event = factor(c("a", "b", "c")), year = c(2007, 2011, 2007),
    category = factor(c("None", "wind", "None")))
  s <- sim_gdr_fixture(site_covs = sites, event_keys = sites["event"], return_joint = TRUE,
    provenance = list(lambdaformula = ~category * year, coefficients = c(intercept = log(100)),
      ecological_scale = 300, reference_center = 2, reference_scale = 0.5))
  expect_s4_class(s$observed, "unmarkedFrameGDR")
  expect_equal(dim(s$joint), c(3, 2, 3))
  expect_equal(apply(s$joint, c(1, 2), sum), s$observed@yDistance)
  expect_equal(apply(s$joint, c(1, 3), sum), s$observed@yRemoval)
  expect_equal(rowSums(s$observed@yDistance), s$truth$realized_observed)
  expect_equal(rowSums(s$observed@yRemoval), s$truth$realized_observed)
  expect_true(all(s$truth$realized_abundance >= s$truth$realized_observed))
  expect_equal(s$truth$expected_abundance, c(100, 200, 0) * pi * 0.075^2)
  expect_equal(s$truth$realized_abundance[3], 0)
  expect_equal(s$truth$expected_observed,
    s$truth$expected_abundance * s$truth$distance_probability * s$truth$removal_probability)
  expect_identical(s$observed@siteCovs$category, sites$category)
  expect_equal(s$observed@siteCovs$year, sites$year)
  expect_false(any(names(s$truth) %in% names(s$observed@siteCovs)))
  expect_identical(s$event_keys, sites["event"])
  expect_no_error(invisible(utils::capture.output(methods::show(s$observed))))
  expect_equal(s$design$provenance$ecological_scale, 300)
  expect_equal(s$design$distance_scale, rep(40, 3))
  expect_identical(s$design$density_units, "birds/km2")
  expect_true(validate_unmarked_data(s$observed, s$event_keys, "event",
    interval_units = s$design$interval_units, effort = rep(4, 3))$valid)
  expect_null(sim_gdr_fixture()$joint)
})

test_that("the no-covariate simulation returns the unchanged native frame contract", {
  s <- sim_gdr_fixture()
  native <- unmarked::unmarkedFrameGDR(yDistance = s$observed@yDistance,
    yRemoval = s$observed@yRemoval, numPrimary = 1L,
    dist.breaks = s$design$distance_breaks, unitsIn = s$design$distance_units,
    period.lengths = s$design$interval_lengths)
  expect_identical(s$observed, native)
})

test_that("half-normal annuli agree with independent integration and distance units", {
  s <- sim_gdr_fixture()
  expected <- vapply(1:2, function(j) stats::integrate(function(r)
    2 * r * exp(-r^2 / (2 * 40^2)) / 75^2,
    lower = c(0, 25)[j], upper = c(25, 75)[j])$value, numeric(1))
  expect_equal(s$probabilities$distance, matrix(rep(expected, each = 3), 3, 2), tolerance = 1e-12)
  km <- sim_gdr_fixture(distance_scale = 0.04, dist_breaks = c(0, 0.025, 0.075), units_in = "km")
  expect_equal(km$truth, s$truth)
  expect_equal(km$observed@yDistance, s$observed@yDistance)
  expect_equal(km$observed@yRemoval, s$observed@yRemoval)
  extreme <- sim_gdr_fixture(distance_scale = c(1e-200, 1e200, 40))
  expect_equal(extreme$truth$distance_probability[1:2], c(0, 1))
  expect_equal(extreme$probabilities$distance[2, ], diff(c(0, 25, 75)^2) / 75^2)
})

test_that("first-removal probabilities respect unequal durations and zero/one limits", {
  s <- sim_gdr_fixture()
  expect_equal(s$probabilities$removal[1, ], c(0.3, 0.7 * (1 - 0.7^2), 0.7^3 * 0.3))
  expect_equal(s$truth$removal_probability, rep(1 - 0.7^4, 3))
  p <- rbind(c(0, 0, 0), c(1, 0.3, 0.7), c(0.2, 0.5, 1))
  s <- sim_gdr_fixture(removal = p)
  expect_equal(s$probabilities$removal[1, ], c(0, 0, 0))
  expect_equal(s$probabilities$removal[2, ], c(1, 0, 0))
  expect_equal(s$probabilities$removal[3, ], c(0.2, 0.8 * 0.75, 0.8 * 0.25))
  expect_equal(s$truth$realized_observed[1], 0)
  one <- sim_gdr_fixture(density = 100, distance_scale = 40, removal = 0.3,
    dist_breaks = c(0, 75), interval_lengths = 1, return_joint = TRUE)
  expect_equal(dim(one$joint), c(1, 1, 1))
  expect_equal(dim(one$observed@yDistance), c(1, 1))
})

test_that("GDR seeded simulation restores existing and absent RNG state", {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (existed) original <- get(".Random.seed", .GlobalEnv)
  on.exit({
    if (existed) assign(".Random.seed", original, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(42); before <- .Random.seed; kind <- RNGkind()
  a <- sim_gdr_fixture(return_joint = TRUE)
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind(), kind)
  expect_identical(sim_gdr_fixture(return_joint = TRUE), a)
  rm(".Random.seed", envir = .GlobalEnv)
  b <- sim_gdr_fixture()
  expect_false(exists(".Random.seed", .GlobalEnv, inherits = FALSE))
  expect_identical(sim_gdr_fixture(), b)
})

test_that("bounded Poisson thinning moments agree with their analytical targets", {
  # Predeclared software distribution check: two seeds, 20,000 independent
  # sites each; six analytic standard errors, with no seed selection or refits.
  n <- 20000L; density <- 8 / (pi * 0.075^2)
  for (seed in c(2281L, 2282L)) {
    s <- sim_gdr_fixture(density = rep(density, n), seed = seed)
    for (nm in c("realized_abundance", "realized_observed")) {
      mu <- if (nm == "realized_abundance") 8 else s$truth$expected_observed[1]
      y <- s$truth[[nm]]
      expect_lt(abs(mean(y) - mu), 6 * sqrt(mu / n))
      expect_lt(abs(stats::var(y) - mu), 6 * sqrt(mu / n + 2 * mu^2 / (n - 1)))
    }
    expect_equal(rowSums(s$observed@yDistance), rowSums(s$observed@yRemoval))
  }
})

test_that("GDR generating probabilities and moments match a fixed native oracle", {
  # One zero-iteration native C fixture supplies exact specified coefficients.
  # It is a probability/simulation oracle, never an inferential or AIC winner.
  n <- 50L; density <- 8 / (pi * 0.075^2)
  s <- sim_gdr_fixture(density = rep(density, n))
  native <- unmarked::gdistremoval(data = s$observed, output = "density",
    unitsOut = "kmsq", keyfun = "halfnorm", mixture = "P", engine = "C",
    K = 80, starts = c(log(density), log(40), qlogis(0.3)),
    se = FALSE, control = list(maxit = 0))
  expect_equal(as.numeric(unmarked::coef(native)), c(log(density), log(40), qlogis(0.3)))
  p <- unmarked::getP(native)
  expect_equal(p$dist[, , 1], s$probabilities$distance, tolerance = 1e-10)
  expect_equal(p$rem[, , 1], s$probabilities$removal, tolerance = 1e-12)
  expect_equal(unmarked::predict(native, "lambda", level = NULL)$Predicted, rep(density, n))
  set.seed(2283L)
  draws <- unmarked::simulate(native, nsim = 200)
  totals <- unlist(lapply(draws, function(x) rowSums(x$yDistance)), use.names = FALSE)
  expected <- s$truth$expected_observed[1]
  expect_lt(abs(mean(totals) - expected), 6 * sqrt(expected / length(totals)))
  expect_true(all(vapply(draws, function(x)
    identical(rowSums(x$yDistance), rowSums(x$yRemoval)), logical(1))))
  # Native zero abundance still uses a proper nonzero conditional detection law.
  native@estimates@estimates$lambda@estimates[] <- -Inf
  zero <- unmarked::simulate(native, nsim = 1)[[1]]
  expect_equal(sum(zero$yDistance), 0)
  expect_equal(sum(zero$yRemoval), 0)
})

test_that("native character-to-factor diagnostics are retained for supplied covariates", {
  expect_warning(s <- sim_gdr_fixture(site_covs = data.frame(category = c("None", "wind", "None"))),
    "contains characters")
  expect_identical(as.character(s$observed@siteCovs$category), c("None", "wind", "None"))
})

test_that("GDR simulator rejects unsupported or ambiguous generating inputs", {
  expect_error(sim_gdr_fixture(density = -1), "density")
  expect_error(sim_gdr_fixture(distance_scale = 0), "distance_scale")
  expect_error(sim_gdr_fixture(distance_scale = 1:2), "one value per site")
  expect_error(sim_gdr_fixture(removal = 1.1), "probabilities")
  expect_error(sim_gdr_fixture(removal = matrix(0.2, 2, 3)), "matrix dimensions")
  expect_error(sim_gdr_fixture(dist_breaks = c(1, 25, 75)), "dist_breaks")
  expect_error(sim_gdr_fixture(interval_lengths = c(1, 1.5)), "integer durations")
  expect_error(sim_gdr_fixture(interval_units = ""), "interval_units")
  expect_error(sim_gdr_fixture(seed = 1.5), "integer seed")
  expect_error(sim_gdr_fixture(site_covs = data.frame(x = 1:2)), "site_covs")
  expect_error(sim_gdr_fixture(return_joint = NA), "return_joint")
  expect_error(sim_gdr_fixture(provenance = list(1)), "named list")
})

test_that("simulation printing works independently of the no-covariate native frame printer", {
  s <- sim_gdr_fixture()
  expect_s3_class(s, "msr_gdr_simulation")
  text <- utils::capture.output(returned <- print(s))
  expect_match(paste(text, collapse = "\n"), "3 survey events")
  expect_match(paste(text, collapse = "\n"), "birds/km2")
  expect_identical(returned, s)
})
