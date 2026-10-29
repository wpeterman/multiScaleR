testthat::test_that("spatial reference inputs match both exact frozen pilots", {
  for (nm in c("null", "smooth")) {
    case <- sr_case(nm)
    old <- jsonlite::read_json(file.path(previous_study, paste0("pilot-", nm, ".json")), simplifyVector = TRUE)
    testthat::expect_identical(case$input_hash, old$input_hash)
    testthat::expect_s3_class(case$object, "msr_spatial_scale_prototype")
    testthat::expect_false(inherits(case$object, "multiScaleR_bayes"))
    testthat::expect_equal(unname(case$truth), old$comparison$truth)
    testthat::expect_identical(case$variables, old$comparison$variable)
  }
})

testthat::test_that("private projection reconstruction matches native predictions and chain order", {
  case <- sr_case("smooth")
  n <- 20L * 4L
  p <- case$object$parameters
  columns <- unlist(lapply(names(p), function(nm) {
    q <- length(p[[nm]])
    if (!q) character() else if (q == 1L) nm else paste0(nm, "[", seq_len(q), "]")
  }))
  flat <- matrix(rep(unlist(p), each = n), nrow = n)
  flat <- flat + .3 * sin(matrix(seq_along(flat), nrow = n))
  raw <- array(flat, c(20L, 4L, ncol(flat)), dimnames = list(NULL, NULL, columns))
  before <- serialize(case, NULL); rng <- .Random.seed
  derived <- sr_derived_array(raw, case)
  unpacked <- sr_unpack(raw, p)
  for (i in c(1L, 20L, 21L, 40L, 61L, 80L)) {
    parameters <- lapply(unpacked, function(x) if (ncol(x)) unname(x[i, ]) else numeric())
    expected <- multiScaleR:::.msr_spatial_scale_predict(case$object, case$data,
      parameters = parameters, scale = "log")
    iteration <- ((i - 1L) %% 20L) + 1L; chain <- ((i - 1L) %/% 20L) + 1L
    testthat::expect_identical(names(expected), rownames(case$data))
    testthat::expect_equal(unname(derived[iteration, chain, 1:72]), unname(expected))
  }
  # A nonidentity Cholesky factor must be included, rather than assuming IID bases.
  block <- case$object$specification$scales$habitat$random$.spatial_field$gaussian_block
  block$chol <- chol(matrix(c(2, .2, .1, .2, 1, .15, .1, .15, .8), 3L))
  block$chol <- t(block$chol)
  case$object$specification$scales$habitat$random$.spatial_field$gaussian_block <- block
  correlated <- sr_derived_array(raw, case)
  parameters <- lapply(unpacked, function(x) if (ncol(x)) unname(x[21L, ]) else numeric())
  expected <- multiScaleR:::.msr_spatial_scale_predict(case$object, case$data,
    parameters = parameters, scale = "log")
  testthat::expect_identical(names(expected), rownames(case$data))
  testthat::expect_equal(unname(correlated[1L, 2L, 1:72]), unname(expected))
  testthat::expect_false(isTRUE(all.equal(correlated, derived)))
  case$object$specification$scales$habitat$random$.spatial_field$gaussian_block$chol <-
    unserialize(before)$object$specification$scales$habitat$random$.spatial_field$gaussian_block$chol
  testthat::expect_identical(serialize(case, NULL), before)
  testthat::expect_identical(.Random.seed, rng)
  threshold <- median(derived[, , 1L])
  indicator <- 1 * (derived[, , 1L] <= threshold)
  event <- sr_event(derived[, , 1L], threshold)
  testthat::expect_equal(event$estimate, mean(indicator))
  testthat::expect_equal(event$mcse, posterior::mcse_mean(indicator))
  testthat::expect_equal(event$ess_mean, posterior::ess_mean(indicator))
  testthat::expect_equal(event$rhat, posterior::rhat(indicator))
  testthat::expect_equal(event$chain_probabilities, unname(colMeans(indicator)))
  quantiles <- sr_quantiles(derived)
  testthat::expect_equal(quantiles$threshold[1:3], unname(quantile(derived[, , 1L], c(.025, .5, .975))))
  testthat::expect_equal(quantiles$quantile_mcse[1:3], unname(posterior::mcse_quantile(derived[, , 1L], c(.025, .5, .975))))
  functional <- sr_functional(derived, quantiles[1:3, ])
  testthat::expect_equal(functional$mcse[2L], event$mcse)
  testthat::expect_equal(unname(as.numeric(functional[2L, paste0("chain", 1:4)])), event$chain_probabilities)
  constant <- sr_event(matrix(0, 20L, 4L), 0)
  testthat::expect_identical(constant$estimate, 1)
  testthat::expect_false(constant$precision_pass)
  testthat::expect_true(is.na(constant$mcse))
  nonfinite <- matrix(0, 20L, 4L); nonfinite[1L, 1L] <- NA
  testthat::expect_match(sr_event(nonfinite, 0)$reason, "no iterations dropped")
  testthat::expect_false(sr_event(nonfinite, 0)$precision_pass)
  testthat::expect_error(sr_unpack(raw[, , -1L, drop = FALSE], p), "parameter layout")
  testthat::expect_error(sr_unpack(matrix(0, 2, 2), p), "dimensions/names")
})
