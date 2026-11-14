#' Glossary of terms used in multiScaleR
#'
#' @description
#' Plain-language definitions of the terms used across the multiScaleR help
#' pages and vignettes, in alphabetical order. Each entry names the page that
#' explains the idea in more depth.
#'
#' @details
#' * **1D and 2D distance**: two ways to turn sigma into a distance. The 2D
#'   distance (the default) is the radius of the circle around a point that
#'   holds a set proportion of an ideal kernel's weight on a plane. The 1D
#'   distance is measured along a line through the point and is shorter. See
#'   [kernel_dist()].
#' * **Backend**: the computational route [multiScale_optim()] uses to search
#'   over scales, chosen with its `backend` argument. See
#'   [multiScaleR-fitting-methods].
#' * **Binned inputs**: distance-binned summaries that [kernel_prep()] stores
#'   when `bin = TRUE` (the default). They let the search compute
#'   kernel-weighted means without visiting every cell, which is much faster.
#' * **Buffer**: the set of raster cells within `max_D` of a survey point. Only
#'   these cells contribute to the point's kernel-weighted covariates.
#' * **Choice set**: in a step-selection model, one used step together with
#'   the available steps it is compared against. See
#'   `vignette("v10-step-selection", package = "multiScaleR")`.
#' * **Conditional uncertainty**: uncertainty that treats the estimated scales
#'   as known. Most coefficient standard errors and prediction intervals are
#'   conditional. See [multiScaleR-prediction].
#' * **Credible interval**: a Bayesian interval that contains a stated
#'   posterior probability (for example 95 percent) for a quantity, given the
#'   model, the data, and the priors. See [multiScale_bayes()].
#' * **Cross-validation**: estimating out-of-sample predictive performance by
#'   fitting the model to part of the data (the training folds) and scoring its
#'   predictions on the rest (the held-out fold). *Spatial folds* hold out
#'   blocks of nearby points, and *grouped folds* hold out whole groups such as
#'   sites or years. See [cross_validate()].
#' * **Distance-removal model**: a model for counts recorded in distance
#'   classes and removal intervals that separates abundance from detection.
#'   multiScaleR uses `unmarked`'s generalized distance-removal (GDR) models.
#'   See `vignette("v12-distance-removal", package = "multiScaleR")`.
#' * **Effective distance**: the distance that contains a set proportion
#'   (by default 90 percent) of a kernel's weight. It converts sigma into a
#'   distance in map units that is easy to report and compare across kernels.
#'   See [kernel_dist()].
#' * **Expected response**: the mean response a model predicts for given
#'   covariate values, such as an expected count or an occurrence probability.
#'   It differs from an *observation replicate*, a simulated new observation
#'   that also includes observation noise. See [multiScaleR-prediction].
#' * **Exponential power kernel** (`"expow"`): a kernel with a shape parameter
#'   in addition to sigma, which lets the weight curve range from sharply
#'   peaked to nearly flat-topped. See [multiScaleR-kernels].
#' * **Fixed kernel**: a kernel that gives every cell within a radius the same
#'   weight and cells beyond it no weight. Its scale is a hard radius.
#' * **Fixed-radius metric**: a landscape or surface metric computed at a
#'   radius you set, so it has no scale to estimate and enters the model as an
#'   ordinary covariate. See [msr_vars()].
#' * **Gaussian kernel**: the default kernel, whose weights decline with
#'   distance following a normal (bell-shaped) curve with standard deviation
#'   sigma. See [multiScaleR-kernels].
#' * **Hard radius**: a scale at which a cell is either fully in or fully out
#'   of the calculation. Landscape metrics, unweighted surface metrics, and the
#'   fixed kernel use hard radii. Their likelihood changes in steps, so they get
#'   no Hessian standard error.
#' * **Hessian**: the matrix of second derivatives of the log-likelihood at its
#'   maximum. Its inverse approximates the variance of the estimates and gives
#'   Wald standard errors.
#' * **Joint uncertainty**: uncertainty that lets the scales and the
#'   coefficients vary together, so that intervals also reflect uncertainty
#'   about the scales. Only some backends provide it. See
#'   [multiScaleR-fitting-methods].
#' * **Kernel**: a function that assigns a weight to each raster cell
#'   according to its distance from a survey point. Kernel-weighted covariates
#'   summarize the landscape around each point. See [multiScaleR-kernels].
#' * **Kernel-weighted mean**: the weighted average of a raster's cell values
#'   around a point, using kernel weights. It is the default covariate type,
#'   set with [kernel_var()].
#' * **Landscape metric**: a summary of the composition or arrangement of
#'   land-cover classes within a radius, such as the proportion of a class or
#'   an edge density, set with [landscape_var()].
#' * **Laplace approximation**: a fast approximation that treats a
#'   distribution (a posterior, or the integral over random effects) as
#'   Gaussian around its peak. See
#'   `vignette("v15-bayesian-computation", package = "multiScaleR")`.
#' * **Lean mode**: preparing inputs with `store_cell_data = FALSE`, which
#'   keeps only the binned summaries to save memory. It works only when every
#'   covariate is a kernel-weighted mean.
#' * **Link scale and response scale**: a generalized linear model works on
#'   the link scale (for example log counts or log odds) and reports
#'   predictions on the response scale (counts or probabilities) after
#'   back-transforming with the inverse link. See [multiScaleR-prediction].
#' * **Marginal likelihood**: the likelihood after integrating out random
#'   quantities, such as random effects or a between-group standard deviation
#'   of scale.
#' * **`max_D` and `min_D`**: the largest and smallest scales the search may
#'   consider, in map units. Set `max_D` in [kernel_prep()] above the largest
#'   plausible scale of effect; `min_D` is set from the raster resolution.
#' * **Negative exponential kernel** (`"exp"`): a kernel whose weights decline
#'   exponentially with distance, giving more weight to nearby cells and a
#'   longer tail than the Gaussian kernel. See [multiScaleR-kernels].
#' * **Posterior**: in a Bayesian model, the probability distribution of the
#'   parameters given the data and the priors. See [multiScale_bayes()].
#' * **Profile likelihood**: the log-likelihood as a function of one parameter,
#'   with every other parameter re-estimated at each value. Profile-likelihood
#'   intervals follow the shape of the likelihood instead of assuming it is
#'   symmetric. See [profile_sigma()].
#' * **Scale of effect**: the distance over which a landscape variable
#'   influences the response. multiScaleR estimates it as the sigma of a
#'   kernel, which [kernel_dist()] converts to an effective distance.
#' * **Shape**: the second parameter of the exponential power kernel, which
#'   controls how sharply the weights fall off.
#' * **Sigma**: the scale parameter of a kernel, in map units. A larger sigma
#'   spreads the weight over a wider area. Sigma is not itself the distance
#'   to report; see *Effective distance*.
#' * **Spatial response field**: a spatially correlated random effect that
#'   captures residual spatial pattern the covariates leave unexplained. See
#'   `vignette("v11-spatial-fields", package = "multiScaleR")`.
#' * **SQ**: the kernel-weighted root-mean-square deviation of a continuous
#'   surface around its local mean, a measure of surface roughness. See
#'   [surface_var()].
#' * **Standardization**: centering a covariate on its mean and dividing by its
#'   standard deviation. multiScaleR standardizes each optimized covariate at
#'   every candidate scale by default, so a coefficient describes the effect of
#'   a one standard deviation change.
#' * **Step-selection model**: a model of animal movement that compares each
#'   used step with available steps through conditional logistic regression.
#'   See `vignette("v10-step-selection", package = "multiScaleR")`.
#' * **Surface metric**: a summary of the texture of a continuous surface
#'   around a point, such as its roughness or skewness, set with
#'   [surface_var()].
#' * **Wald interval**: an interval built from an estimate plus or minus a
#'   multiple of its standard error. It assumes a symmetric, nearly quadratic
#'   likelihood near the estimate, which may not hold for scales.
#'
#' @seealso [multiScaleR-package] for an overview of the package, and
#'   `vignette("v01-quickstart", package = "multiScaleR")` to get started.
#' @name multiScaleR-glossary
NULL
