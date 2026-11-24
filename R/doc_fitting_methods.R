#' Fitting methods for estimating scales of effect
#'
#' @description
#' [multiScale_optim()] estimates scales of effect through one of several
#' *backends*, chosen with its `backend` argument. A backend is the
#' computational route that searches over kernel scales and fits the model at
#' each candidate. Backends differ in the models they accept, the parameters
#' they estimate together, and the uncertainty they can report. This page
#' states the requirements and reporting rules of each backend. To choose a
#' backend for your model, start with
#' `vignette("v07-fitting-methods", package = "multiScaleR")`, which compares
#' them on one data set.
#'
#' @details
#' # Compare the backends
#'
#' | `backend` | Models it fits | Estimated in one search | Uncertainty for the scales |
#' |:----------|:---------------|:------------------------|:---------------------------|
#' | `"refit"` (default) | Any model that `stats::update()`, the `unmarked` equivalent, or a `refit_fn` can refit | Kernel scales and shapes. The model's own parameters are refit at every candidate | Hessian (Wald) standard errors for smooth kernels. Profile likelihood with `summary(fit, profile = TRUE)` |
#' | `"profile_glm"` | Gaussian identity, Poisson log, and binomial logit `glm()` fits, and NB2 log-link fits | Kernel scales and NB2 theta. Coefficients are solved at every candidate | Joint Wald limits when the checks below pass |
#' | `"joint_nb2"` | NB2 log-link fits, including the exact Poisson boundary | Kernel scales, coefficients, and dispersion | Joint Wald limits when the checks pass. None at the Poisson boundary |
#' | `"joint_mixed"` | Poisson and NB2 `glmmTMB` fits with one random intercept | Kernel scales and all model parameters, with the random intercept integrated out | Wald limits for each scale that passes its checks. Profile likelihood with `summary(fit, profile = TRUE)`. Coefficient uncertainty is conditional on the scales |
#' | `"conditional_logit"` | Step-selection models from `survival::clogit()` or `amt::fit_clogit()` | Smooth kernel scales and conditional-logistic coefficients | None. Coefficient uncertainty is conditional on the scales |
#'
#' A *conditional* standard error treats the estimated scales as known. A
#' *joint* standard error lets the scales and the coefficients vary together,
#' so it also carries uncertainty about the scales.
#' [multiScale_bayes()] is a separate route that samples or approximates a
#' posterior instead. See [multiScaleR_bayes-class] and the Bayesian vignettes.
#'
#' # The default refit backend
#'
#' ## How the search works
#'
#' The default backend uses the L-BFGS-B algorithm (a bounded quasi-Newton
#' search) to maximize the log-likelihood over sigma, and over shape for the
#' `"expow"` kernel. At each candidate scale it refits the model, so the
#' regression coefficients take their fitted values for that scale. Candidate
#' scales that leave a covariate undefined at a modeled observation are
#' rejected, so every candidate is compared on the same sample. The final
#' `sample_size` entry of [diagnostics()] reports the prepared and fitted row
#' counts. For `unmarked` data with repeated sites, it also reports the number
#' of distinct prepared spatial sites.
#'
#' ## Standard errors and intervals
#'
#' Standard errors are Hessian-based approximations from this outer search.
#' A nonfinite, singular, or non-positive-definite Hessian gives missing
#' standard errors and a retained diagnostic, and the raw Hessian is kept.
#' So does a Hessian built partly from a failed refit: a refit that fails
#' returns a large penalty, and a finite-difference point that lands on one
#' gives an absurdly large curvature. This case also gives a warning, because
#' it means the model is hard to refit near the estimate.
#'
#' A smooth kernel or shape parameter within `1e-5` of a search bound, on the
#' log scale, invalidates the Wald approximation for that parameter. Its
#' standard error and Wald limits are withheld. The other smooth parameters keep
#' Wald standard errors from their own block of the Hessian, which holds the
#' bounded parameters fixed at their bounds;
#' `diagnostics(fit)$refit_hessian$conditional_on` names those bounded
#' parameters. The raw Hessian, the original curvature status, and a
#' table of parameter bounds stay in `diagnostics(fit)$refit_hessian`. A lower
#' bound of zero is flagged only when it is reached. This flag is numerical. It
#' does not establish that a scale is identified or give a calibrated interval.
#' Native coefficient covariance still conditions on the selected kernels, and
#' likelihoods and parameter counts are unchanged.
#'
#' Hard-radius covariates (landscape metrics, unweighted surface metrics, and
#' every scale under `kernel = "fixed"`) have a likelihood that changes in
#' steps as the radius grows. They receive no Hessian standard error or Wald
#' limit. Their raw Hessian rows remain for compatibility and diagnostics but
#' are not inverted. See [multiScaleR-kernels].
#'
#' Summary methods report profile-likelihood intervals for sigma once
#' [profile_sigma()] has been run on the fit, and Hessian-based intervals
#' otherwise. `summary(fit, profile = TRUE)` computes profile-likelihood
#' intervals directly.
#'
#' # Profile hard radii over a grid
#'
#' The `radius_grid` argument is experimental, so its behavior may change in a
#' later release without a deprecation period. Supply it to the refit backend
#' for a likelihood profile over hard radii. At each combination of radii, the
#' model refits its coefficients and any native nuisance parameters on the same
#' observations. When Gaussian, exponential, or exponential-power kernels supply
#' smooth predictors, a bounded search starts afresh at each combination and
#' refits every smooth bandwidth and applicable shape. Use `par` to give the
#' first smooth start in the usual parameter layout.
#'
#' The largest eligible native log-likelihood selects the combination. Exact
#' ties favor the smallest radii, in model predictor order. Inspect failed
#' combinations and likelihood plateaus in `diagnostics(fit)$radius_profile`.
#' Hard radii have no standard error or confidence interval. Smooth standard
#' errors and Wald limits, when available, condition on the selected radii.
#'
#' The grid refits serially. It accepts ordinary `lm()`, `glm()`, and
#' `MASS::glm.nb()` fits, and complete one-case conditional-logistic models,
#' each with a retained model frame. Fixed predictors, including fixed-radius
#' metrics, can accompany the optimized predictors. Retain raw cells and
#' matching model row IDs, and omit `join_by`, custom refit functions, and start
#' screening. For step selection, supply one case and at least one available
#' alternative per stratum, with an unweighted, unpenalized native likelihood.
#'
#' # Analytic GLM backends
#'
#' `"profile_glm"` and `"joint_nb2"` evaluate the likelihood and its derivatives
#' directly rather than refitting the model through `update()`. This makes
#' joint uncertainty possible.
#'
#' ## Models and data they accept
#'
#' Both need a smooth kernel (Gaussian, exponential, or exponential power) and
#' a `glm` or `MASS::glm.nb` object that retains its model frame. Fit a Gaussian
#' model with `glm(family = gaussian())`, because an `lm()` fit is rejected.
#' `"profile_glm"` supports the Gaussian identity, Poisson log, binomial logit,
#' and NB2 log links. NB2 has variance \eqn{\mu + \mu^2/\theta}, and theta is
#' estimated alongside the kernel parameters even when the starting GLM used a
#' fixed theta. `"joint_nb2"` requires an NB2 log-link model.
#'
#' The retained model-frame rows define the fitting sample. Their row names
#' must match unique prepared point IDs, so preserve them through any
#' subsetting. Kernel means and standard deviations use only these rows, and no
#' observations are dropped during the search.
#'
#' A binomial response may be:
#' * numeric `0`/`1` or logical `FALSE`/`TRUE` (Bernoulli observations);
#' * a two-level factor, with success in its second level;
#' * `cbind(successes, failures)` with positive trial totals. Logical count
#'   columns stand for zero and one counts without changing the trial totals.
#'
#' Response-scale predictions are success probabilities. Positive whole-number
#' frequency weights are allowed for count and binomial models. Gaussian
#' weights must all equal one. Proportion responses and precision weights are
#' not supported.
#'
#' ## Formula terms
#'
#' Numeric kernel predictors may enter native R product, crossing, and nesting
#' terms, including factor interactions and products of optimized predictors.
#' The native terms and fitted contrast coding define the coefficient layout,
#' and every kernel, coefficient, and cross derivative enters the joint
#' uncertainty.
#'
#' Scalar arithmetic on kernel values inside `I()`, and the base functions
#' `exp`, `log`, `sqrt`, `sin`, `cos`, and `tan`, keep their native coefficient
#' basis and analytic derivatives. Their inputs must be optimized kernel values
#' or literal constants, and a nonfinite result is an error.
#'
#' Cubic `splines::ns()` and `splines::bs()` terms also work when you write out
#' literal fixed `knots` and two increasing `Boundary.knots`. All basis columns
#' share the ecological scale of their source, including in factor
#' interactions. Knots are in the units of the modeled predictor. The analytic
#' GLM, conditional-logistic, and joint mixed searches keep their
#' candidate-specific standardization and its derivatives, so fixed
#' standardized knots are not fixed thresholds of raw cover. Bayesian targets
#' keep their frozen reference transformation instead: with a fixed reference,
#' convert raw knots and boundaries with its center and SD before writing them
#' into the formula. Basis values and their derivatives update at each
#' candidate scale, while the knots, the coefficient basis, and any priors stay
#' fixed. Only cubic bases with distinct, strictly interior knots are
#' supported. Adaptive `df` recipes, repeated or boundary knots, and other
#' degrees need the refit backend. Natural-spline tails are linear on the link
#' scale. These are unpenalized bases of fixed complexity, with no
#' smoothing-variance uncertainty and no monotonicity constraint.
#'
#' GAMs keep their smooth terms and penalties only through the refit backend.
#' They are not treated as ordinary unpenalized GLMs here.
#'
#' Not supported: fixed inputs inside a dynamic expression, adaptive dynamic
#' bases, nondifferentiable or arbitrary functions, kernel-dependent offsets,
#' hard-radius metrics, random effects, and zero inflation. Put offsets in the
#' formula with `offset()` rather than in the separate `offset` argument.
#'
#' ## How the search works
#'
#' Both backends run full multi-start searches in serial by default. Use
#' `n_cores = 1`, `start_strategy = "single"`, and `profile_control`. They do
#' not accept `refit_fn` or `join_by`. Every winning start must pass both the
#' optimizer's convergence check and a projected-gradient check.
#'
#' For `"profile_glm"`, regression coefficients are solved by iteratively
#' reweighted least squares (IRLS) at each candidate. A full-rank NB2 fit whose
#' iterations fail receives a coefficient refinement that uses the analytic
#' score and observed information. The Gaussian residual standard deviation is
#' profiled at its maximum likelihood estimate. Independent native likelihoods
#' have verified the fitted objective in simulations of the supported families.
#' This backend is experimental, and the calibration of its scale intervals is
#' not established. Simulations found Wald undercoverage even among fits that
#' passed the curvature checks. For single-kernel NB2 models,
#' [scale_inference()] examines the full likelihood range and can assess
#' individual scale values with a parametric bootstrap.
#'
#' `"joint_nb2"` uses L-BFGS-B to estimate the kernel parameters, the
#' coefficients, and \eqn{\rho = \log(1 + \alpha)} together, where
#' \eqn{\alpha = 1/\theta}. At rho = 0 the likelihood is exactly Poisson.
#' Each candidate must also agree, within `1e-5` negative log-likelihood units,
#' with an independently optimized conditional coefficient fit. Start records
#' keep failures, projected gradients, conditional likelihood gaps, and
#' warnings. `dispersion_est` adds `alpha` and `exact_poisson`, and a boundary
#' estimate has `theta = Inf`. The model is still NB2, so its parameter count
#' includes dispersion at the boundary. Its conditional coefficient summary
#' uses the Poisson GLM summary at that boundary; see `dispersion_est` for the
#' NB2 boundary details. The `profile_fit` parameter vector holds the log kernel
#' parameters followed by rho, and `joint_par` adds the regression
#' coefficients. Interior joint covariance from `vcov()` uses log(theta)
#' coordinates, as for `"profile_glm"`.
#'
#' For NB2 models in either backend, `profile_control$coefficient_start =
#' "warm"` reuses conditional regression coefficients within each search. The
#' default is `"cold"`. Warm continuation resets at every independent start
#' and falls back to cold initialization if its score or rank checks fail.
#' Final likelihood checks always use independently initialized coefficients.
#' Counts and fallback reasons are stored in `profile_fit$coefficient_audit`.
#'
#' ## Run joint_nb2 starts in parallel
#'
#' For `"joint_nb2"` only, `profile_control$worker_pool` spreads the
#' independent starts over a PSOCK cluster you create with
#' `parallel::makePSOCKcluster()`. It does not reduce the number of starts or
#' change any uncertainty rule. Keep `n_cores = 1`, because the pool size sets
#' the number of workers. Workers must load the same installed package code as
#' the main session. Each start keeps its cold terminal audit, and results keep
#' their original start order. Worker coefficient audits are stored in
#' `profile_fit$worker_coefficient_audits`, and `coefficient_audit` then covers
#' the final evaluation in the main session. Later profile and bootstrap
#' inference stays serial.
#'
#' You own the pool. Create it once, pass
#' `profile_control = list(worker_pool = pool)` to each call, and close it with
#' `parallel::stopCluster()` when finished, including after errors (use
#' `on.exit()` inside your own function). Pool connections are not stored in
#' the fitted object. Separate processes duplicate the data and package memory,
#' so check [estimate_multiscale_ram()] first. Pool startup also checks the
#' binned data against 75 percent of physical RAM, allowing four copies of the
#' data and 1 GiB per process for the main session and the workers. This is a
#' planning estimate, not a guarantee of peak memory. Startup can outweigh the
#' gain for short jobs. Set math-library thread limits before creating workers
#' to avoid nested parallel computation.
#'
#' ## Uncertainty rules
#'
#' Joint observed information includes kernel, coefficient, and dispersion
#' cross derivatives. Default scale limits are normal approximations on the
#' log scale. Wald uncertainty is withheld when any of these holds:
#' * a parameter sits at a search bound;
#' * the information is not positive, or its scaled condition number is at
#'   least 1e8;
#' * a log-kernel standard error exceeds one;
#' * competing starts lie more than 0.5 log units apart within a 95 percent
#'   likelihood threshold;
#' * a Poisson fit at the selected scales comes within that same likelihood
#'   threshold of the NB2 fit;
#' * the fit sits on the exact Poisson boundary.
#'
#' These are screening rules, not a proof that the scales are identified. The
#' Poisson comparison is conservative and is not a formal model test;
#' reoptimizing the Poisson scales could reveal more compatibility. A separate
#' near-Poisson diagnostic flags `max(mu/theta) < 0.01`. A finite upper theta
#' bound does not establish a finite dispersion maximum. The marginal
#' covariance and Hessian are missing whenever joint uncertainty is
#' unavailable.
#'
#' A kernel parameter at a search bound does not withhold
#' joint uncertainty when every other rule passes for the remaining
#' parameters. Joint covariance and joint prediction uncertainty then hold the
#' bounded kernel parameters fixed at their bounds. They come from the observed
#' information with those parameters removed, which must be positive and well
#' conditioned, with every remaining log-kernel standard error at most one.
#' Competing modes, a dispersion bound, a Poisson-compatible NB2 fit, or every
#' kernel parameter at a bound still withholds them. `vcov(fit)` returns zero
#' rows and columns for the held parameters, the marginal Hessian is missing
#' for them, and `diagnostics(fit)$profile_glm$joint_conditional_on` names
#' them. `uncertainty_available` stays `FALSE`. The rule is the same for
#' `"profile_glm"` and `"joint_nb2"`.
#'
#' The scale table applies the kernel-specific rules one parameter at a time.
#' A kernel parameter loses its standard error and Wald limits only when it sits
#' at a search bound, its own log-scale standard error exceeds one, or competing
#' starts move it by more than 0.5 log units. Bounded kernel parameters are
#' dropped from the observed information, so the remaining standard errors hold
#' them fixed at their bounds, and `diagnostics(fit)$profile_glm` names them in
#' `scale_se_conditional_on`. Unusable information, a dispersion bound, a
#' Poisson-compatible NB2 fit, or the exact Poisson boundary still withholds
#' every scale standard error.
#'
#' Predictions with `uncertainty = "conditional"` at the Poisson boundary
#' condition on both the kernel parameters and the boundary dispersion. They
#' say nothing about uncertainty in either. See [multiScaleR-prediction].
#'
#' The conditional `opt_mod` keeps the standard GLM methods, but
#' `update(opt_mod)` is disabled because its original call refers to the
#' kernel values from before the search. Refit through [multiScale_optim()].
#'
#' `summary(fit, profile = TRUE)` reoptimizes the other kernel parameters and
#' the NB2 dispersion while it profiles each sigma. It reports the connected
#' likelihood region around the selected solution. Limits that reach a search
#' bound are labeled `"profile-bound-limited"` and are not fully resolved
#' confidence limits. [profile_sigma()] instead holds the other kernel
#' parameters fixed while it reestimates the coefficients and NB2 dispersion.
#' Neither is guaranteed to find likelihood modes the search never visited.
#'
#' ## Compare with an ordinary GLM
#'
#' At the same fixed scales, with the same rows, design, offset, and response
#' coding, the analytic binomial likelihood is the ordinary binomial GLM
#' likelihood. You can check the retained conditional GLM against a separate
#' `glm()` fit on the selected predictors. Do not compare the joint intervals
#' with conditional GLM intervals, because the joint intervals also carry
#' scale uncertainty. For a covariance comparison, use the same response form
#' and weights. A vector response and a two-column response can give slightly
#' different working covariance even when the coefficients and likelihood
#' agree. [aic_tab()] and [bic_tab()] treat logical and numeric binomial
#' responses as the same zero and one observations when the retained rows,
#' trial totals, and likelihood weights agree.
#'
#' # Joint mixed backend
#'
#' `backend = "joint_mixed"` needs the optional packages RTMB, TMB, and
#' glmmTMB. It estimates the kernel scales and the model parameters together,
#' integrating the random intercept out with the Laplace approximation. It fits
#' log-link Poisson and NB2 `glmmTMB` models with one conditional random
#' intercept, constant dispersion, and optional fixed-effects zero inflation.
#' Use Gaussian or exponential kernels, retained model row names that match the
#' prepared point IDs, and `n_cores = 1`. Numeric spatial predictors may enter
#' native R crossing, nesting, and factor interaction terms.
#'
#' Not supported: random slopes, more than one grouping term, spatial terms in
#' the zero-inflation or dispersion formulas, non-unit weights, priors,
#' parameter maps, and REML. Raw formula variables must be recoverable from the
#' fitted model, so put offsets in the formula with `offset()`.
#'
#' Each start uses a fresh native fit, a joint Hessian check, and an
#' independent native refit at its end point whose likelihood must agree within
#' `1e-5`. Failed starts stay in `mixed_fit$start_results`. If every start
#' fails, the error condition carries `start_results`. No Hessian check is
#' bypassed.
#'
#' The result's `opt_mod` is the final native `glmmTMB` fit, and `mixed_fit`
#' records the log scales, bounds, starts, every outcome, and the selected
#' start. Narrow raw prediction data keep fixed-transform inputs for effect
#' curves once the temporary refit data are gone. `diagnostics(fit)$joint_mixed`
#' reports eligibility and scale-boundary flags. As with the analytic GLM
#' backends, an estimate whose distance approaches `max_D` records the
#' `max_distance` diagnostic and warning code, and `fit$next_run` suggests a
#' follow-up `max_D` and starts.
#'
#' Scale standard errors come from the joint Hessian of the selected start,
#' which TMB computes for the Laplace-approximated likelihood with every model
#' parameter free. `fit$scale_est$SE` is sigma times the standard error of
#' log(sigma), and `summary(fit)` reports Wald limits on the log scale, as for
#' `"profile_glm"`. The checks match that backend, one scale at a time. A scale
#' within `1e-5` of a search bound on the log scale loses its standard error,
#' and the other scales hold it fixed at its bound. A scale also loses its
#' standard error when its log-scale standard error exceeds one, or when
#' another start within the 95 percent likelihood threshold moves it by more
#' than 0.5 log units. A Hessian that is not positive definite, or unusable
#' covariance, withholds every scale standard error.
#' `diagnostics(fit)$joint_mixed` gives the status, the reason for each scale,
#' and the bounded scales held fixed. Legacy scale Hessians stay `NA`.
#'
#' Treat these limits with caution. In a simulation with one Gaussian kernel
#' mean and 200 points in 40 sites, 95 percent Wald intervals held the true
#' sigma in 77 percent to 97 percent of interior fits, depending on the scale
#' and response family. They fell short where the likelihood was sharply
#' peaked near the estimate but stayed close to its maximum over a wider range
#' of scales. The
#' limits describe local curvature only, and coefficient and prediction
#' covariance still condition on the selected scales. Use
#' `predict(fit, se.fit = TRUE, uncertainty = "conditional")` or
#' `predict(fit$opt_mod, ...)`. This backend does not extend the
#' single-predictor NB2 contract of [scale_inference()].
#'
#' `summary(fit, profile = TRUE)` reports profile-likelihood limits instead.
#' For each scale in turn, it holds sigma at candidate values and refits the
#' native `glmmTMB` model, the same new fit that checks each start. The
#' coefficients, the random-intercept standard deviation, and any NB2
#' dispersion are re-estimated, and the random intercept is integrated out
#' as in the fit. With several scales, the other scales are re-estimated too,
#' on the joint likelihood from the optimum of each eligible start, and a
#' native refit checks that value. As for the analytic GLM backends, the
#' limits bound the connected region around the estimate where the
#' likelihood-ratio statistic stays below the 95 percent chi-square quantile
#' with one degree of freedom. A limit that reaches a search bound is labeled
#' `"profile-bound-limited"` and is not a fully resolved confidence limit.
#' Every scale gets profile limits, including one whose standard error is
#' withheld, for example at a search bound. Profiling stops with an error when
#' a refit fails or finds a better optimum than the fit, and it needs the
#' optional packages and the refit data stored in `mixed_fit`. In the same
#' simulation, profile intervals held the true sigma in 842 of 900 data sets
#' (89 percent to 96 percent by scale and family), against 821 of 899 for the
#' Wald intervals. Where Wald coverage was lowest, NB2 at a true sigma of 240,
#' the profile intervals held it in 133 of 150 data sets and the Wald
#' intervals in 115 of 149. Profile limits therefore describe scale
#' uncertainty better than the Wald limits, but they can still fall short of
#' nominal coverage, as they did in that scenario. Each profile took a median
#' of 2 to 5 seconds there. [profile_sigma()] instead holds the other scales
#' at their estimates.
#'
#' # Conditional-logistic backend
#'
#' `backend = "conditional_logit"` fits step-selection models with the exact
#' likelihood of a finite choice set: one used step and at least one available
#' step per set. See `vignette("v10-step-selection", package = "multiScaleR")`.
#'
#' ## Models and data it accepts
#'
#' Use a `survival::clogit()` model fitted with `model = TRUE`, or an
#' `amt::fit_clogit()` wrapper, with Gaussian or exponential kernels and one
#' bare `strata(choice_id)` term. A quoted choice-set ID is still one bare
#' variable, and special-call text inside a quoted predictor name is not read
#' as a Cox special term. IDs must identify sets uniquely across animals and
#' bursts, and row names must match the prepared point IDs. Prepare complete
#' sets before fitting. Omitted missing rows, subsets, weights, robust
#' clustering, and special Cox terms are rejected. Keep the available steps
#' fixed throughout the search.
#'
#' Numeric kernel predictors may enter native R interaction, crossing, and
#' nesting terms, including crosses with factors or other optimized
#' predictors. Fixed transformations and offsets keep their native meaning and
#' saved training values. Kernel-only scalar arithmetic and the supported
#' smooth base functions follow the same rules as the analytic GLM backends.
#'
#' Optimized kernel-weighted surface SQ predictors are supported when the
#' matching raw buffers are retained. They use the weighted root-mean-square
#' deviation, not the sample SD within a hard radius. Their scales are
#' estimated separately, and prediction needs the same source and metric
#' specifications and the training normalization. Fixed-radius metrics can
#' enter as fixed predictors. Their original prepared numeric columns must be
#' retained in the native model frame, and new predictions need those columns
#' unless `include_fixed = TRUE` verifies and reuses their original
#' normalization (see [prepare_prediction()]).
#'
#' Not supported: hard-radius landscape and surface metrics with estimated
#' radii, adaptive dynamic bases, arbitrary or nondifferentiable functions,
#' fixed context inside a dynamic expression, and kernel-dependent offsets.
#'
#' ## How the search works
#'
#' The backend runs independent L-BFGS-B starts in serial, with cold
#' coefficient initialization and a cold native `survival` refit at every end
#' point. A start is eligible when it converges, has a small projected score,
#' agrees in likelihood and coefficients, has well-conditioned coefficient
#' information, and agrees with the native covariance. Failed attempts and
#' warnings stay explicit in `profile_fit$start_results`. If every start fails,
#' the error condition carries `start_results`. `profile_fit$starts` stores the
#' log-scale starting values, and `winner` identifies the selected record.
#' After the original joint point passes the native agreement checks, cold
#' native coefficients may refine its terminal score if the objective does not
#' increase and the score improves. Both scores and the refinement flag are
#' kept, and the eligibility threshold is unchanged.
#'
#' ## What the fit keeps
#'
#' The final `opt_mod` is the native conditional fit, with `amt` metadata when
#' supplied. If the original fitting call no longer exposes its raw data,
#' retained model-frame values support fixed transformations, factors, and
#' compound offsets without inverting them. Cold native refits keep the
#' original prediction transformations, including learned polynomial bases.
#' Prediction for new data still needs the raw inputs named in the formula.
#' The internal prediction record omits duplicate fitting data to save space;
#' the original model stays in `fitted_mod_original`. Updating the conditional
#' fit needs explicit data with kernel covariates that are already
#' standardized; otherwise rerun the original multiscale fit.
#'
#' Scale standard errors and joint prediction uncertainty are withheld.
#' Coefficient covariance conditions on the selected scales and does not
#' correct for temporal dependence among steps.
#'
#' ## Movement terms
#'
#' The original `amt` proposal distributions are unchanged. Apply `amt`'s
#' movement-distribution updates to `opt_mod` when they suit the movement
#' terms in your model. Movement terms must match the available-step proposal.
#' The backend does not add them, check a sampling interval, or calibrate
#' population inference.
#'
#' # Covariate types and transformations
#'
#' The analytic backends support optimized kernel means and kernel-weighted
#' surface SQ (the weighted root-mean-square deviation), alone or together.
#' Weighted SQ uses retained raw cell buffers with Gaussian or exponential
#' kernels. Estimating the radius of any other surface metric needs the refit
#' backend. Kernel means need `bin = TRUE` in [kernel_prep()]. Models with only
#' surface metrics need explicit preparation row IDs and raw buffers, without
#' binned summaries. See [msr_vars()] for how each covariate type is defined.
#'
#' A fixed-radius metric enters the native numeric design as an ordinary
#' predictor, with no scale parameter, in every backend and in
#' [multiScale_bayes()]. Its original prepared values and normalization must
#' stay recoverable from the native model frame. If every covariate is fixed,
#' there is no scale to estimate, so fit the model directly to the prepared
#' columns.
#'
#' Under the `"expow"` kernel, a hard-radius predictor has no shape parameter.
#' Its reported shape and shape standard error are `NA`, and model selection
#' leaves the unused shape position out of the parameter count. The starting
#' parameter vector keeps its full layout, so give a finite positive
#' placeholder for that position when you supply `par`.
#'
#' By default (`covariate_transform = "candidate"`), each optimized covariate
#' is centered and scaled again at every candidate scale, so a coefficient
#' describes a change of one standard deviation. The refit backend also accepts
#' `covariate_transform = "raw"`, which keeps optimized kernel and metric
#' values in their own units. The fit records centers of zero and scales of one
#' and reuses them for prediction, effect curves, and spatial
#' cross-validation; they are never re-estimated from a prediction batch. Fit
#' the starting model in those raw units, including any fixed spline knots.
#' Raw mode does not change upstream surface reference normalization or the
#' transforms of fixed metrics, and Bayesian targets are unaffected.
#'
#' The Bayesian route also accepts weighted SQ, with a fixed pooled reference
#' normalization and proper priors. Neither the candidate standardization nor
#' the raw transform of [multiScale_optim()] is substituted into a Bayesian
#' target.
#'
#' # Distance-removal models
#'
#' Native distance-removal (GDR) models from `unmarked` use the refit backend.
#' [make_gdr_refitter()] provides an opt-in native C or TMB refit function
#' when you need independent starts and retained audits. It freezes the
#' evaluated settings and keeps the final native covariance. TMB needs explicit
#' `tmb_settings`, starts in constructor order including log-SDs, and native
#' joint covariance checks. Its intermediate fits still compute `sdreport`.
#'
#' Identified native GDR response fields ([fit_gdr_response_field()]) need a
#' serial scale search: `n_cores = 1` and `PSOCK = FALSE`. A `NULL` core count
#' selects one for this route. Sending native fits to parallel workers has not
#' been qualified.
#'
#' Only active GDR component formulas identify spatial predictors.
#' Availability is inactive for a single primary period, and the distance
#' formula is inactive for a uniform key. Optimized static site covariates must
#' not share names with columns in the native yearly or observation covariate
#' tables. Dot formulas use the original component's available training
#' columns.
#'
#' Native GDR models with active formula offsets or TMB zero inflation are
#' rejected before fitting, model selection, or prediction. In the audited
#' `unmarked` 1.5.2 provider, offsets did not enter the fitted likelihood and
#' the TMB zero-inflated Poisson likelihood ignored its zero-inflation
#' parameter. multiScaleR does not remove offsets, switch engines, change
#' parameter counts, or substitute uncertainty estimates for these cases.
#' Offset-free C zero-inflated Poisson models keep their native likelihood and
#' parameter count. See `vignette("v12-distance-removal", package = "multiScaleR")`.
#'
#' # Control settings
#'
#' `profile_control` is a named list. Entries a backend does not accept are an
#' error.
#'
#' | Entry | Backends | Default | What it controls |
#' |:------|:---------|:--------|:-----------------|
#' | `n_starts` | analytic GLM, `joint_mixed`, `conditional_logit` | `6` | Number of deterministic kernel starts |
#' | `starts` | analytic GLM, `joint_mixed`, `conditional_logit` | `NULL` | Matrix of explicit starts, one row per start (see below) |
#' | `iter.max` | analytic GLM, `conditional_logit` | `500` | Iteration limit for each L-BFGS-B search |
#' | `eval.max` | analytic GLM | `1000` | Evaluation limit for later profile and inference refits, which use conditional fitting |
#' | `gradient_tolerance` | analytic GLM, `conditional_logit` | `0.001` | Largest projected gradient a winning start may have |
#' | `kernel_engine` | analytic GLM, `conditional_logit` | `"reference"` | `"stream"` computes kernel values and derivatives with less copying |
#' | `theta_bounds` | analytic GLM | `c(1e-4, 1e5)`; `c(1e-4, Inf)` for `joint_nb2` | Search bounds for NB2 theta |
#' | `theta_start` | analytic GLM | `NULL` | Starting NB2 theta |
#' | `coefficient_start` | analytic GLM, NB2 only | `"cold"` | `"warm"` reuses coefficients within a search |
#' | `worker_pool` | `joint_nb2` only | `NULL` | PSOCK cluster for running starts in parallel |
#'
#' Each row of `starts` holds positive kernel scales in internal distance units
#' (map distance divided by `unit_conv`), followed by shapes for `"expow"`, in
#' model predictor order. Explicit starts replace the deterministic scale
#' starts, and `par`, if supplied, is added first. All independent starts are
#' kept. NB2 dispersion starts cycle through the initial theta, 0.3, 2, and 20
#' within the bounds. For `"profile_glm"`, theta bounds must be positive,
#' increasing, and at most 1e6. For `"joint_nb2"`, the upper bound must stay
#' infinite and the lower bound must be at least 1e-8. Each `"joint_nb2"`
#' kernel start is tried at the exact Poisson boundary and at an interior
#' dispersion, giving 12 fits with the default six kernel starts.
#'
#' The `"stream"` kernel engine traverses the prepared matrices contiguously,
#' copies them less, and caches response-only NB2 terms. It keeps the same
#' starts, likelihood, and uncertainty checks. It uses extra temporary vectors
#' proportional to the number of observations and can speed up larger binned
#' inputs. For smooth weighted SQ predictors under `"conditional_logit"`, it
#' also computes values and scale derivatives in compiled code, reusing work
#' buffers across sites, with temporary storage proportional to the largest
#' retained site buffer.
#'
#' @seealso [multiScale_optim()] for the arguments, [multiScaleR-class] for the
#'   fitted object, [multiScaleR-kernels] for kernels and hard radii,
#'   [multiScaleR-prediction] for prediction uncertainty, and
#'   `vignette("v07-fitting-methods", package = "multiScaleR")` and
#'   `vignette("v09-estimates-and-intervals", package = "multiScaleR")`.
#' @name multiScaleR-fitting-methods
NULL
