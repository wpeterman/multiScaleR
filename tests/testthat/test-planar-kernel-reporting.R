kernel_plot_layer <- function(plot, geom) {
  which(vapply(plot$layers, function(layer) inherits(layer$geom, geom), logical(1)))
}

test_that("default and legacy plot markers agree with all kernel definitions", {
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    beta <- if (kernel == "expow") 1.5 else NULL
    for (dimension in c("2d", "1d")) {
      p <- plot_kernel(sigma = 300, kernel = kernel, beta = beta,
                        dimension = dimension)
      line <- kernel_plot_layer(p, "GeomVline")
      expected <- k_dist(300, prob = 0.9, kernel = kernel, beta = beta,
                          dimension = dimension)
      expect_equal(p$layers[[line]]$data$xintercept, expected)
      label <- p$labels$caption
      expect_match(label, if (dimension == "2d") "ideal planar weight \\(2D\\)" else
                     "line-based weight \\(1D\\)")
      expect_match(label, "Sigma: 300")
      expect_gte(max(p$data$dist), expected)
      expect_equal(min(p$data$dist), 0)
    }
    default <- plot_kernel(sigma = 300, kernel = kernel, beta = beta)
    expect_equal(default$layers[[kernel_plot_layer(default, "GeomVline")]]$data$xintercept,
                 k_dist(300, prob = 0.9, kernel = kernel, beta = beta, dimension = "2d"))
  }
})

test_that("small scales and high-proportion markers are not rounded away or clipped", {
  for (dimension in c("2d", "1d")) {
    for (kernel in c("gaussian", "exp", "expow")) {
      beta <- if (kernel == "expow") 0.5 else NULL
      p <- plot_kernel(sigma = 0.01, kernel = kernel, beta = beta,
                        prob = 0.99999, dimension = dimension)
      marker <- p$layers[[kernel_plot_layer(p, "GeomVline")]]$data$xintercept
      expect_gt(marker, 0)
      expect_gte(max(p$data$dist), marker)
      expect_warning(ggplot2::ggplot_build(p), NA)
    }
  }
  # This exact quantile remains defined beyond the old 10-sigma root bracket.
  r <- k_dist(1, prob = 0.99999, kernel = "expow", beta = 0.5, dimension = "1d")
  expect_gt(r, 10)
  expect_equal(pgamma(sqrt(r), shape = 2), 0.99999)
})

test_that("exact exponential-power distances integrate to the requested weight", {
  for (beta in c(0.5, 1, 1.5, 2, 5)) {
    for (dimension in c("2d", "1d")) {
      radius <- k_dist(3, prob = 0.9, kernel = "expow", beta = beta,
                        dimension = dimension)
      # Normalize analytically, then integrate independently of the quantile.
      power <- if (dimension == "2d") 1 else 0
      norm <- 3^(power + 1) * gamma((power + 1) / beta) / beta
      weight <- integrate(function(d) d^power * exp(-(d / 3)^beta),
                           0, radius, rel.tol = 1e-9)$value / norm
      expect_equal(weight, 0.9, tolerance = 1e-7)
    }
  }
})

test_that("summary and printing select 2D or legacy 1D without modifying the fit", {
  fit <- make_core_fixture()$opt
  before <- fit
  planar <- summary(fit)
  legacy <- summary(fit, dimension = "1d", prob = 0.95)
  expect_identical(planar$opt_dist, planar$opt_dist_2d)
  expect_identical(legacy$opt_dist, legacy$opt_dist_1d)
  expect_identical(legacy$opt_dist, legacy$opt_distance)
  expect_equal(legacy$opt_dist, kernel_dist(fit, prob = 0.95, dimension = "1d"))
  expect_identical(fit, before)
  for (dimension in c("2d", "1d")) {
    expected <- if (dimension == "2d") "2D, radius" else "1D, along"
    printed <- capture.output(print(planar, dimension = dimension))
    expect_match(printed[grep("selected", printed)], expected)
    direct <- capture.output(print(fit, dimension = dimension, prob = 0.95))
    expect_match(direct[grep("selected", direct)], expected)
    expect_true(any(grepl("95%", direct, fixed = TRUE)))
  }
  old <- planar
  old$opt_dist <- old$opt_dist_1d
  old$opt_dist_1d <- old$dimension <- NULL
  printed <- capture.output(print(old))
  expect_match(printed[grep("selected", printed)], "1D, along")
})

test_that("fitted plot markers and intervals use the chosen dimension", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  fit <- make_core_fixture()$opt
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    fit$kernel_inputs$kernel <- kernel
    fit$shape_est <- if (kernel == "expow") data.frame(Mean = 1.5, SE = 0.1,
                                                      row.names = "cont1") else NULL
    for (dimension in c("2d", "1d")) {
      tab <- kernel_dist(fit, prob = 0.95, dimension = dimension)
      p <- plot(fit, prob = 0.95, dimension = dimension)[[1]]
      line <- p$layers[[kernel_plot_layer(p, "GeomVline")]]$data$xintercept
      rect <- ggplot2::ggplot_build(p)$data[[kernel_plot_layer(p, "GeomRect")]]
      expect_equal(line, tab$Mean)
      expect_equal(unique(rect$xmin), tab[, 2])
      expect_equal(unique(rect$xmax), tab[, 3])
      expect_gte(max(p$data$dist), tab[, 3])
      expect_match(p$labels$caption,
                   if (dimension == "2d") "ideal planar weight" else "line-based weight")
    }
  }
})

test_that("hard-radius plots show buffers without invented weight claims or intervals", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  fit <- make_core_fixture()$opt
  fit$kernel_inputs$scale_vars$type[] <- "landscape"
  fit$scale_est$SE <- NA_real_
  for (dimension in c("2d", "1d")) {
    p <- plot(fit, prob = 0.8, dimension = dimension)[[1]]
    expect_length(kernel_plot_layer(p, "GeomRect"), 0)
    expect_equal(p$layers[[kernel_plot_layer(p, "GeomVline")]]$data$xintercept,
                 round(fit$scale_est$Mean, 2))
    label <- p$labels$caption
    expect_match(label, "Buffer radius")
    expect_false(grepl("80%", label, fixed = TRUE))
    expect_true(all(p$data$wt[p$data$dist >= fit$scale_est$Mean] == 0))
    expect_length(unique(p$data$wt[p$data$dist < fit$scale_est$Mean]), 1)
    expect_warning(ggplot2::ggplot_build(p), NA)
  }
})

test_that("fitted plots omit missing weighted scales without leaving empty list entries", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  fit <- make_core_fixture()$opt
  fit$scale_est <- rbind(missing = data.frame(Mean = 10, SE = NaN), fit$scale_est)
  plots <- plot(fit)
  expect_length(plots, 1)
  expect_s3_class(plots[[1]], "ggplot")
  expect_identical(plots[[1]]$labels$title, "cont1")
  fit$scale_est$SE[] <- NaN
  expect_identical(plot(fit), list())
})

test_that("plot switches suppress layers and validate probability and dimension", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  fit <- make_core_fixture()$opt
  plain <- plot(fit, scale_dist = FALSE)[[1]]
  expect_length(kernel_plot_layer(plain, "GeomVline"), 0)
  expect_length(kernel_plot_layer(plain, "GeomRect"), 0)
  expect_length(kernel_plot_layer(plain, "GeomText"), 0)
  expect_null(plain$labels$caption)
  expect_length(kernel_plot_layer(plot(fit, add_label = FALSE)[[1]], "GeomText"), 0)
  expect_null(plot(fit, add_label = FALSE)[[1]]$labels$caption)
  expect_length(kernel_plot_layer(plot_kernel(sigma = 10, kernel = "gaussian",
                                              scale_dist = FALSE), "GeomVline"), 0)
  for (prob in list(0, 1, NA_real_, Inf, c(0.8, 0.9), "90%")) {
    expect_error(plot(fit, prob = prob), "prob")
    expect_error(plot_kernel(sigma = 10, kernel = "gaussian", prob = prob), "prob")
  }
  expect_error(plot(fit, dimension = "3d"), "arg")
  expect_error(summary(fit, dimension = "3d"), "arg")
  expect_error(print(fit, dimension = "3d"), "arg")
  expect_error(plot_kernel(sigma = 10, kernel = "gaussian", dimension = "3d"), "arg")
})

test_that("computational cutoffs retain explicit 1D semantics", {
  fix <- make_core_fixture()
  check_distance <- function(model, prob = 0.9, dimension, ...) {
    expect_identical(dimension, "1d")
    if (missing(model)) kernel_dist(prob = prob, dimension = dimension, ...) else
      kernel_dist(model, prob = prob, dimension = dimension, ...)
  }
  # Keep the captured function outside the mocked namespace to avoid recursion.
  original <- kernel_dist
  environment(check_distance) <- list2env(list(kernel_dist = original), parent = environment())
  with_mocked_bindings({
    sim_dat(alpha = 0.5, beta = 0.75, sigma = 40, n_points = 5,
             raster_stack = fix$rs, user_seed = 42)
    sim_dat_unmarked(alpha = 0.5, beta = 0.75, sigma = 40, n_points = 5,
                      n_surv = 3, raster_stack = fix$rs, user_seed = 42)
    kernel_scale.raster(fix$rs, sigma = 40, verbose = FALSE)
    .surface_kernel_weight_matrix(fix$rs, 40, NULL, "gaussian", 0.975)
  }, kernel_dist = check_distance, .package = "multiScaleR")
  expect_equal(diagnostics(fix$opt)$max_distance$effective_distance,
               kernel_dist(fix$opt, dimension = "1d")$Mean, ignore_attr = TRUE)
})
