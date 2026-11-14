#' Bayesian landscape-scale fits
#'
#' @description
#' [multiScale_bayes()] returns an object of class `multiScaleR_bayes`. It is a
#' list that holds the frozen model, the parameter layout, and, once a method
#' has been run, the posterior draws and the record of that computation. This
#' page lists every element, shows where each model parameter lives and on
#' which scale, explains how to convert it to a distance or a multiplier, and
#' names the functions that read the object. The Bayesian interface is
#' experimental. Printing any of these objects reports "Experimental Bayesian
#' landscape-scale model".
#'
#' @details
#' # Elements in every object
#'
#' * `specification`: the frozen model. Rows, landscape bins, formulas,
#'   priors, and the reference standardization are fixed here once and reused
#'   by every method. See "Read the specification".
#' * `parameters`: a named list of seven numeric vectors that defines the
#'   parameter layout and holds starting values. It contains no estimates.
#'   See "Find each parameter".
#' * `call`: the matched call to [multiScale_bayes()].
#' * `inference`: a label for the computation, one of `"prepared"`,
#'   `"full Bayesian NUTS"`, `"full Bayesian slice sampling"`,
#'   `"joint Laplace Gaussian approximation"`,
#'   `"marginal Laplace with scale-SD quadrature"`, or
#'   `"importance-corrected marginal Laplace"`.
#' * `inference_metadata`: a description of what the object estimates. See
#'   "Read the inference metadata". Fits saved by older versions of the
#'   package can lack it, and they remain usable.
#'
#' # Elements added by each method
#'
#' | `method` | Elements added |
#' |:---------|:---------------|
#' | `"prepare"` | None. There are no draws, so the object cannot be used for posterior summaries or prediction |
#' | `"sample"` (NUTS) | `draws`, `sampler`, `sampling`, `diagnostics` |
#' | `"slice"` | `draws`, `chains`, `sampling`, `diagnostics` |
#' | `"laplace"` | `draws`, `laplace`, `diagnostics` |
#' | `"marginal"` | `draws`, `marginal`, `diagnostics` |
#' | `"marginal"` with `correction = "importance"` | `draws`, `marginal`, `importance`, `diagnostics` |
#'
#' * `draws`: a named list of seven numeric matrices, one per parameter block
#'   in the order of `parameters`. Each matrix has one row per retained joint
#'   posterior vector and one column per parameter in the block, on the
#'   internal scales given in "Find each parameter". Row `i` of every block
#'   belongs to the same joint vector, so compute derived quantities from
#'   whole rows. NUTS rows come in permuted order. Slice rows follow the
#'   chains, one chain after another. Laplace and marginal rows are
#'   independent draws from the approximation. Importance-corrected rows are
#'   resampled, with replacement, from the weighted candidates.
#' * `sampler` (NUTS): the `stanfit` object returned by tmbstan.
#' * `chains` (slice): a numeric array of post-warmup iterations by chains by
#'   parameters, with parameter names such as `"beta[1]"` and
#'   `"scale_beta[1]"`. Chain diagnostics and bayesplot chain plots read it.
#' * `sampling` (NUTS and slice): the `chains`, `iter`, `warmup`, `seed`, and
#'   effective `cores` used. NUTS also records `control` and `laplace = FALSE`.
#' * `laplace`, `marginal`, `importance`: records of the approximations. See
#'   "Read the approximation records".
#' * `diagnostics`: the numerical checks of the computation. See "Read the
#'   stored diagnostics".
#'
#' # Find each parameter
#'
#' The seven parameter blocks appear in the same order in `parameters`,
#' `specification$parameters`, `draws`, and every concatenated vector or
#' matrix (the Laplace `mode`, `covariance`, and `hessian`, the importance
#' `candidates`, and the columns of `chains`).
#'
#' | Block | One entry per | What it holds | Scale stored | Interpretable form |
#' |:------|:--------------|:--------------|:-------------|:-------------------|
#' | `beta` | name in `specification$beta_names` | Response coefficients | Link scale | Read directly, or `exp()` for count and binomial models. Summarized in `summary(fit)$coefficients` |
#' | `scale_beta` | column of each scale design | Scale-formula coefficients | Log of distance divided by `specification$unit_conv` for intercepts and `~0 + year` coefficients. Log-distance contrasts for the rest | Add `log(unit_conv)` to the intercepts to get log map units, then `exp()` for a distance. `exp()` of a contrast is a ratio of distances. Summarized in `summary(fit)$scale_coefficients` |
#' | `scale_z` | level of each scale grouping term | Standardized group deviations with a standard-normal prior | None | Multiply by the group SD to get the deviation in log scale |
#' | `log_scale_sd` | scale random-intercept term | Log of the between-group SD of log scale (tau) | Log | `exp()` gives tau. `exp(tau)` is a multiplier on the scale. Summarized in `summary(fit)$group_sd$log_scale_group_sd` |
#' | `response_z` | level of each response grouping term | Standardized random intercepts with a standard-normal prior | None | Multiply by the response group SD to get the intercept on the link scale |
#' | `log_response_sd` | response random-intercept term | Log of the random-intercept SD | Log of a link-scale SD | `exp()`. Summarized in `summary(fit)$group_sd$response_group_sd` |
#' | `log_dispersion` | NB2 or Gaussian model (empty otherwise) | Log NB2 size theta, or log Gaussian residual SD | Log | `exp()`. Summarized in `summary(fit)$dispersion` |
#'
#' The link scale is the log expected count for Poisson and NB2 models, the
#' log odds for binomial models, and response units for Gaussian models.
#'
#' Each landscape predictor's scale for observation `i`, in map units, is
#'
#' \deqn{\log D_i = x_i^\top b + \sum_k \exp(s_k) z_{k, g_k(i)} + \log(u)}{log D[i] = x[i]' b + sum over k of exp(s[k]) z[k, g[k](i)] + log(u)}
#'
#' where \eqn{x_i} is row `i` of the scale design
#' `specification$scales[[name]]$X`, \eqn{b} is
#' `scale_beta[specification$scales[[name]]$beta_index]`, \eqn{s_k} is the
#' entry of `log_scale_sd` for random term \eqn{k}, \eqn{z_{k, g}} is the
#' entry of `scale_z` for group \eqn{g}, and \eqn{u} is
#' `specification$unit_conv`. Within `specification$scales[[name]]$random`,
#' each term records `z_index` (positions in `scale_z`, one per level in
#' `levels`), `index` (each row's level), and `sd_index` (position in
#' `log_scale_sd`). Response random terms in `specification$response$random`
#' use the same fields for `response_z` and `log_response_sd`. This
#' noncentered form, a standard-normal deviation times a group SD, is the
#' coordinate system the samplers and approximations work in.
#'
#' For a shared scale, `exp(scale_beta + log(unit_conv))` is the scale itself.
#' For `~1 + (1 | year)`, the same expression gives the typical scale, the
#' median of the population distribution of annual scales, and tau is
#' `exp(log_scale_sd)`.
#' A year one SD above typical has a scale `exp(tau)` times the typical one.
#' [posterior_scales()] performs these conversions for every distinct
#' combination of scale covariates and groups. Prefer it to rebuilding them
#' by hand.
#'
#' # Read the specification
#'
#' `specification` is a list. Most users need only a few elements.
#'
#' * `covariates`: names of the landscape predictors whose scales are
#'   estimated. `kernel` is `"gaussian"` or `"exp"`.
#' * `unit_conv`: the distance normalizer of the prepared inputs. Internal
#'   log scales are in units of distance divided by `unit_conv`.
#' * `support`: `min_D` and `max_D` of the prepared inputs, in map units. They
#'   are a reference range for [posterior_scales()], not bounds on the prior.
#' * `center`, `spread`: the fixed reference standardization of each
#'   landscape predictor, computed once across fitting rows at the prior
#'   median scale. A response coefficient refers to one `spread` unit of the
#'   kernel-weighted predictor at that scale.
#' * `priors`: the complete prior settings, with defaults filled in.
#'   `scale_median` is a named vector, one entry per landscape predictor.
#' * `beta_names`: the names of the response coefficients, in `beta` order.
#' * `scales`: one list per landscape predictor, holding the scale `formula`,
#'   its design matrix `X`, `beta_index`, `distance_location` (`TRUE` for
#'   coefficients that are log distances rather than contrasts), the prior
#'   mean and SD of each coefficient (`prior_mean`, `prior_sd`, on the
#'   internal log scale), `random` terms, and `labels` with the
#'   `representatives` rows used to report each distinct combination.
#' * `response`: the response formula, its model frame, terms, `random`
#'   terms, and the levels and contrasts retained for prediction.
#' * `family`, `y`, `trials`, `offset`: the response family, the response
#'   vector (successes for binomial data), the binomial trial totals (1 for
#'   other families), and the offset.
#' * `data`, `rows`, `observation_ids`, `row_id`: the observation data, the
#'   prepared landscape row of each observation, its point ID, and the name of
#'   the row-ID column (`NULL` when row names were used). Older saved fits can
#'   lack `row_id`.
#' * `V`, `C`, `delta`: binned value sums, finite-cell counts, and shifted
#'   distance terms for kernel means.
#' * `X_fixed`: the design columns that do not involve a landscape predictor.
#'
#' Some models add `response_design` (interactions and expressions that
#' involve a landscape predictor), `surface` (weighted SQ buffers, their
#' specifications, and the `memory` planning record), or `fixed_reference`
#' (the preparation normalization of declared fixed-radius metrics).
#'
#' # Read the inference metadata
#'
#' `inference_metadata` describes the target and the treatment of uncertainty.
#' It does not report whether a fit passed its checks.
#'
#' * `schema_version`: version of this structure.
#' * `target`: always `"posterior"` for these objects.
#' * `method`: `"prepare"`, `"nuts"`, `"slice"`, `"joint_laplace"`,
#'   `"marginal_laplace"`, or `"importance_corrected"`.
#' * `objective`: `"negative_log_posterior"`. `population_priors` is
#'   `"proper"`, because every parameter has a prior. `response_normalizers`
#'   is `"retained"`.
#' * `integration`: how scale effects, response effects, and population
#'   parameters (group SDs) are treated, such as `"joint_sampling"`,
#'   `"joint_gaussian"`, or `"joint_importance_sampling"`, or
#'   `"not_computed"` for a prepared object.
#' * `solver`: `"none"`, `"nuts"`, `"slice"`, `"optimization"`, or
#'   `"importance_sampling"`.
#' * `uncertainty`: `intervals` is `"marginal_credible"` once fitted and
#'   `"not_provided"` when prepared. `population_parameters` and
#'   `annual_scales` state that group-SD uncertainty is propagated through the
#'   joint posterior. `calibration` is always `"not_established"`, because no
#'   method has demonstrated interval coverage for general data sets.
#' * `diagnostics_field`: the name of the element that holds the checks,
#'   `"diagnostics"`. `metadata_establishes_eligibility` is always `FALSE`.
#'
#' # Read the stored diagnostics
#'
#' Every fit has `diagnostics$passed` (`TRUE`, `FALSE`, or `NA`) and, in most
#' cases, `criteria`, a sentence stating the rule applied. The other fields
#' depend on the method.
#'
#' | Method | Fields | `passed` requires |
#' |:-------|:-------|:------------------|
#' | NUTS | `parameters` (a posterior summary table with `rhat`, `ess_bulk`, `ess_tail`), `divergences`, `max_treedepth_hits`, `ebfmi` (one per chain) | At least 2 chains, Rhat at most 1.01, bulk and tail effective sample size at least 100 per chain, no divergences or tree-depth hits, and E-BFMI at least 0.3 |
#' | Slice | `parameters` | The same chain rules, without the NUTS-only checks. Without the posterior package, `passed` is `NA`, `parameters` is `NULL`, and `reason` explains why |
#' | Laplace | `type`, `gradient_max`, `min_eigenvalue`, `condition_number`, `competing_modes`, `successful_starts`, `starts` | No distinct mode within 2 log-posterior units. A failed optimizer or an invalid Hessian stops the fit instead |
#' | Marginal | `type`, `grid_converged`, `grid_nodes`, `refinement`, `competing_mass`, `upper_tail_mass`, `gradient_max`, `condition_number` | A converged grid, `competing_mass` at most 0.01, and `upper_tail_mass` at most 0.001 |
#' | Importance | `pareto_k`, `ess`, `min_tail_ess`, `max_weight`, `split_cdf_error`, `type`, `candidates`, `min_candidates`, `proposal_df`, `production_candidates`, `proposal_passed` | Pareto k at most 0.7, `ess` at least 1000, `min_tail_ess` at least 400, `max_weight` at most 0.01, and `split_cdf_error` at most 0.05 |
#'
#' In the importance checks, `ess` measures how concentrated the weights
#' are, `min_tail_ess` the precision of tail indicators, `pareto_k` the
#' presence of extreme ratios, `max_weight` dominance by one candidate, and
#' `split_cdf_error` the disagreement between two independent halves of the
#' production candidates. `min_candidates` is the requested minimum.
#'
#' For the samplers, these checks cover the raw parameters. Use
#' [diagnose_bayes_scales()] to check the scales and group SDs you report. For
#' the approximations, a passed check says the arithmetic finished cleanly. It
#' does not say the approximation is close to the posterior.
#' `vignette("v15-bayesian-computation", package = "multiScaleR")` explains
#' each check and what to do when one fails.
#'
#' # Read the approximation records
#'
#' Modes, covariances, Hessians, and candidates in these records use the
#' unconstrained, noncentered coordinates of "Find each parameter", with the
#' blocks concatenated in order.
#'
#' * `laplace`: `mode`, the joint `covariance` and `hessian` at the mode,
#'   every optimizer attempt in `starts`, the `selected_start`, and the number
#'   of `draws`, `seed`, `cores`, and a `coordinates` note. The covariance
#'   keeps correlations among coefficients, scales, and group SDs.
#' * `marginal`: the integration grid over the one scale-group SD. `x` holds
#'   the nodes in the coordinates of the SD's half-normal prior CDF (0 to 1).
#'   `weights`, `density`, `cumulative`, and `log_normalizer` describe the
#'   integrated SD distribution at those nodes. `nodes` holds the conditional
#'   approximation at each node, with its `mode`, `covariance`, `hessian`,
#'   `log_density`, and optimizer `attempts`. The record also keeps
#'   `prior_sd`, the grid `refinement` history, `converged`, `hyper_index`
#'   (the position of the SD), `parameter_indices` (positions of node
#'   parameters in the full layout), and the `draws`, `seed`, `cores`,
#'   `control`, and `coordinates` used. With one response grouping term
#'   integrated inside an importance-corrected fit, the nodes exclude its
#'   latent site effects. The log normalizer describes the numerical
#'   integration. It is not a model-selection criterion.
#' * `importance`: all retained production `candidates` (one row each) with
#'   normalized `weights` and `log_weights` (log posterior minus log proposal
#'   density), and the `resample_ids` into the candidates that built `draws`.
#'   `history` has one row per production stage, with its diagnostics,
#'   whether it `passed`, and whether it met the requested minimum
#'   (`budget_met`). The record also keeps the `seed`, a text description of
#'   the `proposal`, the `adapted_proposal` mixture (its `components` and
#'   their `weights`, where each learned component counts in
#'   `skipped_nodes` the nodes that kept their original geometry),
#'   `pilot_candidates`, `proposal_indices`, the
#'   `conditional_response` flag, `production_candidates` (including
#'   discarded stages), `response_proposal` (fallback count and largest
#'   remaining Newton step), and `proposal_diagnostics`, the original
#'   marginal checks. `target_evaluation` counts every posterior density
#'   evaluation (`evaluations`), invalid automatic-differentiation values
#'   (`ad_invalid`), ordinary rechecks that recovered a finite or zero
#'   density (`recovered_finite`, `recovered_zero`), and values left
#'   `unresolved`, including pilots and discarded stages. Every candidate is
#'   kept, and unresolved densities still fail the importance checks. With
#'   `conditional_beta = TRUE` the record also holds `conditional_beta` and
#'   `beta_proposal`, which gives the starting coefficients (`start`), the
#'   number of `fallbacks` to the coefficients' normal prior, the largest
#'   remaining Newton step (`max_step`), and the `latent_prior_fraction` of
#'   the defensive component.
#'
#' # Functions that read a fit
#'
#' | Function | What it returns |
#' |:---------|:----------------|
#' | `print()` | Family, observations, scale formulas, and the stored diagnostic status |
#' | `summary()` | Coefficients, scales, scale coefficients, group SDs, dispersion, and diagnostics. See [multiScale_bayes()] |
#' | [posterior_scales()] | Each scale in map units or log map units, with mass below `min_D` and above `max_D` |
#' | [plot.multiScaleR_bayes()] | Interval plots of scales or group SDs |
#' | [diagnose_bayes_scales()] | Chain diagnostics and Monte Carlo error for scales and group SDs (needs `chains` or `sampler`) |
#' | [diagnose_bayes_probabilities()] | Probabilities that scales fall below thresholds, with Monte Carlo error |
#' | [diagnose_bayes_laplace()] | A reweighted check of scale probabilities from an uncorrected approximation |
#' | [predict.multiScaleR_bayes()] | Expected responses or linear predictors at fitted or new observations |
#' | [predict_scales()] | Scales at new scale covariates and groups, including unsurveyed years |
#' | [simulate.multiScaleR_bayes()] | Posterior predictive replicate responses |
#' | [simulate_bayes_prior()] | Prior predictive responses and scales (a prepared object is enough) |
#' | [plot_marginal_effects()] | Coefficient effect curves at the fixed reference standardization |
#' | [cross_validate_bayes()] | Refits the model in each fold and scores held-out observations |
#'
#' [aic_tab()] and [bic_tab()] do not accept these objects.
#'
#' @seealso [multiScale_bayes()], [posterior_scales()],
#'   [predict.multiScaleR_bayes()], [diagnose_bayes_scales()];
#'   `vignette("v13-bayesian-scales", package = "multiScaleR")` for fitting and
#'   reading a model,
#'   `vignette("v14-bayesian-prediction", package = "multiScaleR")` for
#'   prediction, and
#'   `vignette("v15-bayesian-computation", package = "multiScaleR")` for the
#'   computation methods and their diagnostics.
#'
#' @examples
#' # A stored slice-sampling fit with one shared scale per landscape layer
#' cache <- readRDS(system.file("extdata", "cache-bayesian-scales.rds",
#'                              package = "multiScaleR"))
#' fit <- cache$fit
#' fit$call                       # the call that produced this fit
#' lengths(fit$parameters)        # 4 coefficients and 2 scale intercepts
#' fit$specification$beta_names   # order of the beta columns
#' dim(fit$chains)                # iterations x chains x parameters
#'
#' # Rebuild sigma for land1, in meters, from whole rows of the draws
#' sm <- fit$specification$scales$land1
#' log_sigma <- fit$draws$scale_beta[, sm$beta_index] +
#'   log(fit$specification$unit_conv)
#' quantile(exp(log_sigma), c(0.025, 0.5, 0.975))
#' posterior_scales(fit)[1, c("median", "lower", "upper")]   # the same values
#'
#' # A stored NUTS fit in which the scale varies among years
#' annual <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'                               package = "multiScaleR"))$fits$pooled
#' year_term <- annual$specification$scales$habitat$random$year
#' tau <- exp(annual$draws$log_scale_sd[, year_term$sd_index])
#' median(tau)        # between-year SD of log sigma
#' median(exp(tau))   # a year one SD above typical has this multiple of sigma
#'
#' @name multiScaleR_bayes-class
NULL
