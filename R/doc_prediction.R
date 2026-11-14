#' Prediction quantities and their uncertainty
#'
#' @description
#' A prediction from a fitted multiscale model is one specific quantity, such
#' as an expected count, a linear predictor, or a detection probability. Its
#' standard error or interval carries one specific kind of uncertainty. This
#' page lists what each prediction route in multiScaleR returns and what its
#' uncertainty includes, so that you can choose the call that answers your
#' question and describe the result correctly. It covers fits from
#' [multiScale_optim()]. For fits from [multiScale_bayes()], see
#' [predict.multiScaleR_bayes()].
#'
#' `vignette("v09-estimates-and-intervals", package = "multiScaleR")` works
#' through these distinctions with examples, and
#' `vignette("v03-map-predictions", package = "multiScaleR")` shows how to map
#' predictions and limit extrapolation.
#'
#' @details
#' # Choose a prediction route
#'
#' Two terms recur below. *Conditional* uncertainty treats the estimated kernel
#' scales (sigma) and shapes as known and carries only the uncertainty in the
#' regression coefficients (and an estimated dispersion, where the model has
#' one). *Joint* uncertainty lets the kernel parameters and the coefficients
#' vary together.
#'
#' | Route | Call | Quantity returned | Standard errors | Intervals |
#' |:------|:-----|:------------------|:----------------|:----------|
#' | Analytic GLM backends (`"profile_glm"`, `"joint_nb2"`) at points | `predict(fit, newdata, kernel_inputs)` | Expected response (`type = "response"`) or linear predictor (`"link"`) | Joint (the default) or conditional | `interval = "confidence"` gives Wald limits for the expected response |
#' | Analytic GLM backends on a raster | `predict(fit, raster)` | Expected response or linear predictor in each cell | Conditional only | None |
#' | Conditional logistic (`"conditional_logit"`) | `predict(fit, newdata, kernel_inputs, type = "link")` | Linear predictor (`"link"`), relative selection score (`"risk"`), or probability within each complete choice set (`"response"`) | Conditional only | None |
#' | Joint mixed model (`"joint_mixed"`) | `predict(fit, newdata, kernel_inputs)` | Native glmmTMB response or link prediction | Conditional only | None |
#' | Default refit backend, most model classes | `predict(fit, newdata, kernel_inputs)` | The provider's prediction: response scale by default, linear predictor for Cox models | Conditional, when the provider supplies them | None from the wrapper. Native `lm` intervals through `native_args` |
#' | Default refit backend, `unmarked` models | `predict(fit, newdata, kernel_inputs, join_by = keys, type = "state")` | One named native component, such as abundance, occupancy, or detection | Conditional | Native confidence limits in the returned table |
#' | Fixed-row snapshot | `predict(prepare_prediction(fit, ...))` | Same as the route it was prepared from | Same | Same |
#' | Effect curves | `plot_marginal_effects(fit)` | Expected response across one covariate, other inputs held at reference values | Conditional, for every backend | 95% Wald band |
#' | Native model | `predict(fit$opt_mod, newdata)` or `terra::predict(raster, fit$opt_mod)` | The provider's own prediction, from covariates you have already standardized | Conditional | The provider's own |
#' | Parameter covariance | `vcov(fit)` | Covariance of the estimated parameters | Joint or conditional (analytic backends) | Not applicable |
#' | Bayesian fit | `predict(bfit, ...)`, `simulate(bfit, ...)` | Posterior expected response, or replicated observations | Full posterior | Credible intervals |
#'
#' Rasters for `predict(fit, raster)` come from [kernel_scale.raster()] with
#' `scale_center = TRUE` and `match_buffer = TRUE`. [predict.multiScaleR()]
#' documents the arguments for the first six routes.
#'
#' # Read the scale of a prediction
#'
#' The *linear predictor*, or link scale, is the intercept plus each
#' coefficient times its covariate, plus any formula offset. The *response
#' scale* applies the model's inverse link function to it. `exp()` inverts a
#' log link and turns a linear predictor into an expected count. `plogis()`
#' inverts a logit link and turns it into a probability. Choose with
#' `type = "link"` or `type = "response"`.
#'
#' To back-transform a link-scale result yourself, apply the inverse link to
#' the estimate and to each confidence limit. Do not apply it to a standard
#' error. A response-scale standard error comes from the delta method,
#' `abs(d mu / d eta) x SE(link)`, where `mu` is the expected response and
#' `eta` is the linear predictor. For a log link this is `mu x SE(link)`.
#' The analytic backends return it with `type = "response", se.fit = TRUE`.
#' Their confidence limits are computed on the link scale and then
#' transformed, so on the response scale they are asymmetric and stay inside
#' the valid range (positive counts, probabilities between 0 and 1).
#'
#' Some models have their own scales:
#' * **Cox and step-selection models.** The linear predictor is a log
#'   relative risk, and `exp()` of it is a relative risk or relative selection
#'   score. Neither is a probability of occurrence or survival. Native Cox
#'   refits use `type = "lp"` (which `type = "link"` also selects) or
#'   `type = "risk"`. The conditional logistic backend uses `"link"` and
#'   `"risk"`, and its `"response"` is a probability normalized within each
#'   supplied choice set.
#' * **Linear models.** For an `lm` refit the link and response scales are
#'   the same, so `type = "link"` returns the response.
#' * **`unmarked` models.** Native prediction back-transforms each component
#'   by default, so abundance is a count and occupancy and detection are
#'   probabilities. Bind `native_args = list(backTransform = FALSE)` in
#'   [prepare_prediction()] to keep the component's link scale.
#'
#' Landscape covariates enter the model standardized by default, so a
#' coefficient describes the change in the linear predictor for one standard
#' deviation of the kernel-weighted covariate. Every prediction route reuses
#' the mean and standard deviation from the fitted data. Never standardize a
#' prediction data set on its own.
#'
#' # Expected response or a new observation
#'
#' The *expected response* is the model's mean at a site: the expected count,
#' or the probability of occurrence. An observation recorded at that site
#' scatters around this mean. That scatter is *observation noise*. Every
#' standard error and interval that [predict.multiScaleR()],
#' [prepare_prediction()], and [plot_marginal_effects()] return describes the
#' expected response or a model component, and none adds observation noise.
#' Two native controls are the exceptions:
#' * For an `lm` refit, `native_args = list(interval = "prediction")` returns
#'   limits for a future observation, including the residual variance.
#' * For a spaMM refit, the native `residVar` and `respVar` variance
#'   attributes are on the response scale and include observation noise.
#'
#' multiScaleR has no `simulate()` method for [multiScale_optim()] fits. For
#' a GLM, the native `simulate(fit$opt_mod)` draws new responses at the fitted
#' rows from the fitted means, holding the scales and coefficients at their
#' estimates. For a Bayesian fit, [simulate.multiScaleR_bayes()] draws
#' replicated observations that carry both posterior uncertainty and
#' observation noise.
#'
#' # Choose conditional or joint uncertainty
#'
#' | Fitting route | Joint uncertainty | Conditional uncertainty |
#' |:--------------|:------------------|:------------------------|
#' | `backend = "profile_glm"` or `"joint_nb2"` | Yes, from `predict()` at points and from `vcov()`, when `fit$uncertainty$available` is `TRUE` | Yes, at points and on rasters |
#' | `backend = "conditional_logit"` | No | Yes |
#' | `backend = "joint_mixed"` | No | Yes, from native glmmTMB |
#' | Default `backend = "refit"` | No | Yes, from the native provider |
#' | [multiScale_bayes()] | The posterior already covers scales and coefficients together | Not applicable |
#'
#' `predict()` defaults to `uncertainty = "joint"`. Every route except the
#' analytic backends at points needs `uncertainty = "conditional"` before it
#' returns a standard error, and stops otherwise rather than silently
#' returning a conditional one. For the analytic backends, joint uncertainty
#' is refused when the fit sits at a parameter bound or fails its
#' identifiability checks. `fit$uncertainty$reasons` says why. Point
#' predictions remain available in that case. At the exact Poisson boundary
#' of a `"joint_nb2"` fit, conditional standard errors and intervals also hold
#' the dispersion at that boundary.
#'
#' Effect curves from [plot_marginal_effects()] and predictions from the
#' native model in `fit$opt_mod` are conditional for every backend, including
#' the analytic ones. A conditional interval answers the narrower question of
#' how uncertain a prediction is if the scales of effect are exactly their
#' estimates. Say so when you report one.
#'
#' # Confidence intervals and prediction intervals
#'
#' A *confidence interval* describes uncertainty about the expected response.
#' A *prediction interval* describes where a single new observation may fall.
#' It adds observation noise and is wider.
#' * The analytic backends return Wald confidence limits with
#'   `interval = "confidence"`, at the confidence `level` you set.
#' * Every other route requires `interval = "none"` in `predict()`. Native
#'   `unmarked` tables carry `lower` and `upper` confidence limits for the
#'   component at the native level, 0.95 unless you bind another `level` in
#'   `native_args`. An `lm` refit returns native mean confidence limits with
#'   `native_args = list(interval = "confidence")` and prediction limits with
#'   `interval = "prediction"`, both with `uncertainty = "conditional"`.
#' * Effect-curve bands from [plot_marginal_effects()] are 95% Wald limits for
#'   the expected response.
#' * Bayesian credible intervals from [predict.multiScaleR_bayes()] describe
#'   the expected response. Intervals from [simulate.multiScaleR_bayes()]
#'   replicates describe new observations.
#'
#' Wald intervals are approximations, and the package does not establish
#' their coverage. `vignette("v09-estimates-and-intervals", package =
#' "multiScaleR")` summarizes how well they performed in simulations.
#'
#' # Random effects and new groups
#'
#' Native mixed models predict with the fitted random effects of known groups
#' by default (glmmTMB and lme4). Setting the random effects to zero gives a
#' reference prediction for a group whose deviation is exactly zero:
#' * glmmTMB, lme4, and spaMM: `re.form = NA` (or `re.form = ~0`).
#' * nlme: `native_args = list(level = 0)`.
#' * `unmarked`: any non-`NULL` `re.form` sets every random effect to zero. It
#'   does not select a subset of effects.
#'
#' A zero-effect prediction is not an average over the random-effect
#' distribution. Under a nonlinear link, the average over groups generally
#' differs from the prediction at zero. None of the frequentist routes integrates over
#' groups or samples effects for new groups. A new group level needs the
#' provider's own control, such as `allow.new.levels = TRUE` for glmmTMB, and
#' then follows the provider's rules. [plot_marginal_effects()] always sets
#' random effects to zero. Bayesian predictions can sample new-group effects,
#' as [predict.multiScaleR_bayes()] describes.
#'
#' Pass `re.form` and other native controls through `...` in
#' [predict.multiScaleR()], or bind them with `native_args` in
#' [prepare_prediction()]. The `effects` field of the [prediction_contract]
#' records the policy that was applied.
#'
#' # Components of unmarked fits
#'
#' An `unmarked` model has several regression components, and a prediction
#' describes one of them. Name it with `type` in [predict.multiScaleR()] and
#' [plot_marginal_effects()], or with `component` in [prepare_prediction()].
#' Point predictions from `occu()` and `pcount()` fits default to `"state"`.
#' Snapshots and raster predictions require an explicit component.
#'
#' | Model | `type` | What it predicts | Rows it uses |
#' |:------|:-------|:-----------------|:-------------|
#' | `unmarked::occu()` | `"state"` | Probability that a site is occupied | Sites |
#' | `unmarked::occu()` | `"det"` | Probability of detection on one survey, given that the site is occupied | Site by survey |
#' | `unmarked::pcount()` | `"state"` | Expected abundance at a site. With `mixture = "ZIP"`, the overall mean, including structural zeros | Sites |
#' | `unmarked::pcount()` | `"det"` | Probability of detecting each individual on one survey | Site by survey |
#' | `unmarked::gdistremoval()` | `"lambda"` | Abundance or density parameter | Sites |
#' | `unmarked::gdistremoval()` | `"phi"` | Availability, active only with more than one primary period | Site by primary period |
#' | `unmarked::gdistremoval()` | `"dist"` | Distance-detection parameter of the fitted key, inactive for a uniform key | Site by primary period |
#' | `unmarked::gdistremoval()` | `"rem"` | Removal parameter | Site by primary period by removal interval |
#'
#' Abundance from `pcount()` counts the animals present, not the animals you
#' would expect to count. For `plot_marginal_effects()`, `type = "lambda"` is
#' an alias for `"state"` in `occu()` and `pcount()` fits, and `"state"` is an
#' alias for `"lambda"` in `gdistremoval()` fits. Other classes keep their
#' native component names, such as `"psi"`, `"col"`, and `"ext"` for dynamic
#' occupancy (initial occupancy, colonization, and extinction). These are
#' distinct from smoothed or projected occupancy trajectories.
#'
#' Native standard errors and confidence limits for every component condition
#' on the selected scales and add no observation noise. A component
#' prediction does not change the fitted model's parameter count, likelihood,
#' or diagnostics. `vignette("v02-user-guide", package = "multiScaleR")`
#' fits and maps N-mixture and occupancy models, and
#' `vignette("v12-distance-removal", package = "multiScaleR")` covers
#' `gdistremoval()`.
#'
#' # Bayesian fits
#'
#' A fit from [multiScale_bayes()] carries posterior draws of the scales and
#' coefficients together, so its predictions need no conditional or joint
#' choice. [predict.multiScaleR_bayes()] returns posterior summaries of the
#' expected response, [simulate.multiScaleR_bayes()] returns replicated
#' observations, and [predict_scales()] returns posterior scales for observed
#' or new groups. `vignette("v14-bayesian-prediction", package =
#' "multiScaleR")` predicts new sites, new years, and totals.
#'
#' # Check what a result describes
#'
#' Snapshots from [prepare_prediction()], native refit predictions, and
#' Bayesian summaries carry a [prediction_contract] that names the quantity,
#' its scale, the random-effect policy, and the uncertainty it includes and
#' omits. Read it when you are unsure which interval you are holding.
#'
#' @seealso [predict.multiScaleR()], [vcov.multiScaleR()],
#'   [prepare_prediction()], [prepare_unmarked_data()],
#'   [plot_marginal_effects()], [prediction_contract],
#'   [predict.multiScaleR_bayes()], [multiScaleR-class][multiScaleR-class],
#'   [multiScaleR-fitting-methods]
#' @name multiScaleR-prediction
NULL
