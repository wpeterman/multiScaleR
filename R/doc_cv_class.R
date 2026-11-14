#' Cross-validation results
#'
#' The cross-validation functions return a list of class `multiScaleR_cv`.
#' This page describes every element of that list, what each score measures
#' and on what scale, and which functions read the result. Several functions
#' add a subclass that changes or adds a few elements, and those differences
#' are described at the end.
#'
#' @details
#' # Which functions return it
#'
#' | Function | Class of the result |
#' |:---|:---|
#' | [cross_validate()], [cross_validate_glm()] (Gaussian, Poisson, binomial) | `multiScaleR_cv` |
#' | [cross_validate()], [cross_validate_glm()] (`"nbinom2"`, `"Gamma"`) | `multiScaleR_distribution_cv`, `multiScaleR_cv` |
#' | [cross_validate_clogit()] | `multiScaleR_choice_cv`, `multiScaleR_cv` |
#' | [cross_validate_field()] | `multiScaleR_field_cv`, `multiScaleR_cv` |
#' | [cross_validate_mixed()] | `multiScaleR_mixed_cv`, `multiScaleR_cv` |
#' | [cross_validate_bayes()] | `multiScaleR_bayes_cv`, `multiScaleR_cv` |
#' | [select_cv()] | the class of its candidates |
#'
#' A *fold* is one split of the data into *analysis* rows, which fit the
#' model, and *assessment* rows, which are held out and scored. An
#' *assessment occurrence* is one held-out row in one fold. A row assessed in
#' two folds has two occurrences.
#'
#' # Read these elements first
#'
#' * `complete`: `TRUE` when every planned fold was fitted and scored. When
#'   any fold fails it is `FALSE`, and every pooled score in `summary` is
#'   `NA`. Never rank a model on the folds that happened to finish.
#' * `summary`: data frame of pooled scores, one row per metric.
#'   * `metric` names the score, and `estimate` is the pooled score (lower
#'     is better).
#'   * `unit` says what one score refers to, `"observation"` or, for the
#'     Brier score, `"trial"`.
#'   * `n_assessment` is the number of planned assessment occurrences, and
#'     `n_scored` the number that entered the score.
#'   * `weight_sum` is the total frequency weight behind the score
#'     (frequency weight x trials for the Brier score).
#'   * `complete` repeats the completion flag.
#' * `metrics`: the same columns plus `fold`, with one row per fold and
#'   metric. A fold that finished keeps its scores even when another fold
#'   failed. The pooled RMSE is not the mean of the fold RMSEs.
#' * `fold_results`: data frame with `fold`, `n_analysis` and `n_assessment`
#'   (rows used to fit and to score), and `status`, either `"scored"` or
#'   `"failed"`.
#' * `errors`: list with one element per fold. It is `NULL` for a scored
#'   fold. For a failed fold it holds `stage`, where the fold stopped (for
#'   example `"initial_glm"`, `"scale_fit"`, `"prediction"`, or
#'   `"training_selection"` in a [select_cv()] result), and `message`, the
#'   reason. A failed scale fit can also carry `start_results` or
#'   `radius_profile`. Failed folds are never retried.
#' * `warnings`: list with one character vector per fold, holding the warning
#'   messages captured while that fold was fitted and predicted. They are
#'   captured rather than printed, so read them here.
#'
#' # Read the scores
#'
#' Every score is a loss, so lower is better. Each assessment occurrence gets
#' its own loss, stored in a column of `predictions` named after the metric.
#' The pooled `estimate` is the frequency-weighted mean of those losses over
#' the eligible occurrences. An occurrence is eligible when its frequency
#' weight is positive and, for binomial data, it has at least one trial.
#' Ineligible occurrences stay in `predictions` with `eligible = FALSE` and
#' contribute nothing.
#'
#' | Metric | Loss for one occurrence | Pooled score and its scale |
#' |:---|:---|:---|
#' | `rmse` | Squared error, (observed - predicted)^2 | Square root of the mean squared error, in response units (for example, animals per survey) |
#' | `mae` | Absolute error | Mean absolute error, in response units |
#' | `neg_log_score` | Negative natural log of the probability (counts, binomial successes) or density (Gamma) that the fitted model gives the observed value | Mean per observation, in natural-log units |
#' | `deviance` | Poisson deviance contribution (Poisson only) | Mean per observation |
#' | `brier` | Squared error of the success probability, averaged over the row's trials (binomial only) | Mean per trial, between 0 and 1 |
#'
#' ## How to read the scores
#'
#' * RMSE and MAE judge only the predicted mean. RMSE weighs large misses
#'   heavily, and MAE treats every miss in proportion to its size. Note that
#'   the `rmse` column of `predictions` holds squared errors, not their roots.
#' * The negative log score judges the whole predicted distribution, so it
#'   also penalizes a model that is too confident. It keeps every constant of
#'   the probability function, such as the binomial coefficient. A
#'   difference of D between two models means the better one gave the
#'   observed outcomes exp(D) times the probability, as a geometric mean per
#'   observation. A difference of 0.1 is a ratio of about 1.11. For Gamma
#'   responses the score is a density, so it depends on the response units
#'   and can be negative.
#' * An observed outcome given zero probability produces an infinite log
#'   score. The value is kept, not clipped.
#' * For NB2 and Gamma, the log score uses theta or the dispersion estimated
#'   from the analysis rows and frozen for scoring.
#' * A score means little alone. Compare candidates fitted to the same folds
#'   with [cv_tab()].
#'
#' # Other elements
#'
#' * `predictions`: data frame with one row per assessment occurrence,
#'   including the occurrences of failed folds, whose predictions and losses
#'   are `NA`.
#'   * `fold` and `row_id` identify the occurrence.
#'   * `observed` is the response (the number of successes for binomial
#'     data), and `trials` is 1 except for `cbind(successes, failures)`
#'     responses.
#'   * `prediction` is the expected response on the response scale (the mean
#'     count, the mean, or the probability of success per trial).
#'   * `weight` is the frequency weight, and `eligible` marks the
#'     occurrences that enter the scores.
#'   * One loss column per metric is named after the metric.
#'   * NB2 and Gamma results add `distribution_family`, the frozen training
#'     `theta` (NB2) or `dispersion` (Gamma), with `NA` where it does not
#'     apply, and `poisson_limit`, which is `TRUE` only when a joint NB2 fit
#'     verified the exact Poisson limit.
#' * `training`: list with one element per fold, describing the fits to the
#'   analysis rows. An element can be empty or partial when its fold failed
#'   early.
#'   * `coefficients`: final regression coefficients on the link scale (log
#'     for Poisson, NB2, and Gamma, logit for binomial, identity for
#'     Gaussian). Use `exp()` on a log-link coefficient for a multiplicative
#'     change in the expected response. Spatial covariates are standardized
#'     with the training mean and SD, so their coefficients refer to a one
#'     training SD change, unless `covariate_transform = "raw"` was used.
#'   * `initial_coefficients`: coefficients of the initial GLM, fitted at the
#'     starting scales.
#'   * `scale_est`: data frame with `Mean` and `SE`, the estimated sigma and
#'     its standard error in map units, one row per optimized covariate.
#'     Similar values across folds suggest that the data identify the scales
#'     well.
#'   * `shape_est`: estimated kernel shape (unitless) for kernels that have
#'     one, otherwise `NULL`.
#'   * `initial_normalization`, `fitted_normalization`: lists of `mean` and
#'     `sd`, the analysis-row centering and scaling of each spatial covariate
#'     at the starting and the fitted scales. Assessment rows reuse them.
#'   * `bin_reference`: the frozen distance bins (`breaks`, representative
#'     distances `bin_rep`, the extraction envelope `domain`, and
#'     `cell_size`). Distances are in units of `max_D`, so 1 equals `max_D`.
#'   * `optimizer`: `par` (the optimized scale parameters on the internal
#'     scale, sigma divided by `max_D`), `value` (the minimized objective),
#'     `convergence` (0 means success), and `message`.
#'   * `start_results`: the table of optimizer starts from
#'     [multiScale_optim()], when available.
#'   * `initial_distribution` (NB2 and Gamma): the nuisance parameter of the
#'     initial fit, in the same form as `distribution`.
#'
#'   [cross_validate_glm()] estimates no scales, so its `training` holds
#'   only `coefficients` (plus the NB2 and Gamma items).
#' * `distribution` (NB2 and Gamma only): list with one element per fold,
#'   holding the frozen training distribution. An element is `NULL` when the
#'   training fit failed.
#'   * `family`, with `theta` (the NB2 size, `Inf` only for a verified
#'     Poisson limit) or `dispersion` (the Gamma dispersion).
#'   * `poisson_limit`, and `estimated`, which is `FALSE` when the parameter
#'     was fixed in advance.
#'   * `convention`, how the value was obtained, for example
#'     `"native_gamma_deviance_scale"` (the training deviance divided by the
#'     sum of the weights).
#'   * `rank`, `native_K`, `logLik`, `nobs`, `training_deviance`, and
#'     `sum_weights` of the training fit.
#' * `diagnostics`: list with one element per fold. `initial_glm` and
#'   `fitted_glm` hold `converged`, `boundary`, and `iter` from the initial
#'   and final GLM (NB2 and Gamma results also record `th.warn` and
#'   `SE.theta` where the native fit has them).
#'   `scale` holds the scale-fit diagnostics described in [diagnostics()].
#'   Check them before ranking models, because a scored fold is not
#'   necessarily a trustworthy one.
#' * `prediction_contracts`: list of the prediction contracts attached by
#'   the prediction method of each fold (see [prediction_contract]). Each
#'   states what the prediction conditions on and which uncertainty it
#'   leaves out. It can be empty when the prediction method attaches none.
#' * `folds`: the split plan as supplied.
#' * `assessment_counts`: named integer vector, one value per input row,
#'   giving how many times the row was assessed (0 for rows never held out).
#' * `targets`: data frame with `row_id`, `observed`, `trials`, and `weight`
#'   for every row used in any split, sorted by row ID. [cv_tab()] and
#'   [select_cv()] use it to confirm that candidates were scored on the same
#'   data. It never enters a fit.
#' * `selection_records`: `NULL` unless `keep_selection = TRUE`. Otherwise a
#'   list with one compact record per fold, read by [select_cv()]. Each
#'   record holds `version`, `controls`, `diagnostics`, metadata `warnings`,
#'   `data`, and an `error` message if the record could not be built. Leave
#'   the records unchanged. They cannot rebuild the fitted models. A
#'   successful `data` holds these items.
#'   * `logLik`, the training log-likelihood (natural log).
#'   * `native_K`, the parameter count of the fitted regression as its
#'     log-likelihood reports it, including an estimated dispersion or theta.
#'   * `ecological_K`, the number of estimated scales and shapes, and `K`,
#'     their total.
#'   * `nobs` and `observation_ids`, the rows that entered the likelihood.
#'   * `targets`, the analysis-row responses, trials, and weights, and
#'     `signature`, a summary of the likelihood used to check that
#'     candidates are comparable.
#'   * `formula` (as text), `family`, and `link`. NB2 and Gamma records add
#'     `nuisance`.
#' * `fits`: `NULL` unless `keep_fits = TRUE`. Otherwise the full model of
#'   each fold, fitted to that fold's analysis rows only.
#' * `raster_groups`: present only for inputs from
#'   [kernel_prep_by_group()]. For each fold, the map labels of the analysis
#'   and assessment rows. These are fixed inputs, not fitted random effects.
#' * `settings`: list with `family`, `link`, `backend` (`"glm"` for
#'   [cross_validate_glm()] and `"training_selection"` for [select_cv()]),
#'   `metrics`, the control lists used, and labels that state what the
#'   scores mean. `normalization = "training_only"` and
#'   `parameter_uncertainty = "not_integrated"`, for example, say that
#'   standardization came from the analysis rows and that the predictions
#'   treat fitted parameters as known. [cv_tab()] checks these labels
#'   before comparing candidates.
#' * `call`: the matched call.
#'
#' A [select_cv()] result also holds `selection`, with the `criterion` and
#' `ties` rule, the `candidates`, the `metrics` kept and `metrics_dropped`,
#' `chosen` (the candidate selected in each fold, `NA` when selection
#' failed), `tables` (the training information-criterion table of each
#' fold), `candidate_settings`, `candidate_warnings`, and
#' `candidate_diagnostics`. Its `settings$controls_from` names the candidate
#' whose controls fill the shared `settings`. NB2 and Gamma selections add
#' `settings$fold_settings`, the settings of the model chosen in each fold.
#'
#' # Functions that read it
#'
#' * `print()` reports how many folds were scored and how many input rows
#'   were assessed, warns when the result is incomplete, and prints
#'   `summary`.
#' * [cv_tab()] ranks complete results on one score, after checking that
#'   every candidate used identical folds, responses, trials, and weights.
#' * [select_cv()] scores a rule that picks a model in each fold by its
#'   training AIC, AICc, or BIC. It needs `selection_records` or `fits`.
#' * [aic_tab()] and [bic_tab()] refuse a cross-validation result, because
#'   held-out scores are not a maximized likelihood.
#'
#' # Conditional choice results (multiScaleR_choice_cv)
#'
#' [cross_validate_clogit()] scores the probability of the used alternative
#' within each held-out choice set, and the elements change accordingly.
#'
#' * `predictions` has one row per alternative. `observed` is 1 for the used
#'   alternative and 0 for the available ones, `choice` is the set ID, and
#'   `trials` and `weight` are 1. `linear_predictor` is uncentered, so only
#'   differences within a set matter. `prediction` is the probability of the
#'   alternative within its set (the probabilities of a set sum to 1), and
#'   `log_probability` is its natural log. Only the used alternative is
#'   `eligible`, and its `neg_log_score` is the negative log probability of
#'   the choice. Log scores are computed directly on the log scale, so they
#'   stay finite when a printed probability rounds to zero.
#' * `summary` and `metrics` hold `neg_log_score` only, averaged over choice
#'   sets, with `unit = "choice_set"`. `n_assessment`, `n_scored`, and
#'   `weight_sum` count sets, and `n_alternatives` counts the alternative
#'   rows. A model without covariates or offsets gives each of k
#'   alternatives a probability of 1/k and scores log(k), for example 1.099
#'   for three alternatives.
#' * `fold_results` counts alternative rows, and `targets` adds `choice`.
#' * `diagnostics` holds `initial` and `native` (iteration count `iter` and
#'   partial-likelihood `method`) and `scale`.
#' * `selection_records` are version 2. Their `nobs` counts choice sets (one
#'   event each), not alternatives, and they record the `method`.
#' * `settings` gives `family = "conditional_choice"`, `link = "softmax"`,
#'   the `method`, and `backend = "clogit"` when no scales were estimated.
#' * Its `print()` method reports the folds scored and prints `summary`.
#'
#' # Spatial field results (multiScaleR_field_cv)
#'
#' [cross_validate_field()] scores each held-out block jointly, so that
#' correlation among neighboring held-out observations counts.
#'
#' * `predictions` holds `fold`, `row_id`, `observed`, `prediction` (the
#'   predicted response mean), `squared_error`, `absolute_error`, and
#'   `eligible`, plus four variance columns.
#'   `new_location_field_variance` is the field variance that would remain
#'   even if the field were known exactly at the training locations, and
#'   `fitted_training_field_variance` is the extra variance that arises
#'   because the field there is itself estimated. Both are on the
#'   linear-predictor scale and, for non-Gaussian families, describe a local
#'   Gaussian approximation. `observation_noise_variance` (response scale)
#'   is the observation variance, averaged over the field for non-Gaussian
#'   families. `observation_predictive_SE` (response scale) is the standard
#'   deviation of a new observation, combining the field and the noise. It is
#'   a predictive SD, not the standard error of the mean.
#' * `joint_scores` holds, for each fold, the block size `n` and the block's
#'   `neg_log_score`. Gaussian blocks add the `quadratic` loss and
#'   `log_determinant`. Laplace blocks add the analysis and combined log
#'   densities, their latent-mode diagnostics, and the gap from the native
#'   analysis likelihood.
#' * `summary` has columns `metric`, `estimate`, `n_assessment`, `unit`, and
#'   `complete`. The log score is the sum of the block scores divided by the
#'   number of assessment observations, with unit
#'   `"joint_block_density_per_observation"` or
#'   `"joint_block_mass_per_observation"` (counts and 0/1 data). RMSE and MAE
#'   have unit `"response"`. There is no `metrics` table, so use
#'   `cv_tab(by_fold = TRUE)` for fold scores.
#' * `training` holds the analysis `row_ids`, `parameters` (`field_SD`, the
#'   correlation `range` in coordinate units, `observation_SD` or the NB2 or
#'   Gamma dispersion, and parameter counts), `coefficients` (link scale),
#'   and, with kernel inputs, `normalization`, `scale_est`, `bin_reference`,
#'   and the number of `refits`.
#' * `selection_records` are always kept. Each holds `logLik`, `K`,
#'   `native_K`, `ecological_K`, `nobs`, `training_ids`, and hashes of the
#'   training responses and coordinates.
#' * `diagnostics` holds the native fit record of each fold, including the
#'   covariance starts tried and the path of its audit file.
#' * `targets` holds `row_id` and `observed`, `coordinates` holds `row_id`,
#'   `x`, and `y`, and `integrity_sha256` is a hash of the recorded results.
#'   [cv_tab()] and [select_cv()] refuse a result whose records were changed
#'   after the run.
#'
#' # Mixed model results (multiScaleR_mixed_cv)
#'
#' [cross_validate_mixed()] scores all held-out rows of a group jointly,
#' because they share one unknown group effect.
#'
#' * `predictions` adds `group`, and `prediction` is the response mean
#'   integrated over the group effect (the success probability per trial for
#'   binomial data).
#' * `joint_scores` holds, for each scored fold, a data frame with `group`,
#'   the number of held-out rows `n`, and the group's joint
#'   `neg_log_score`.
#' * In `summary` and `metrics`, the log score is the sum of the group scores
#'   divided by the number of assessment rows (unit
#'   `"joint_group_density_per_row"`). RMSE and MAE have unit `"response"`,
#'   or `"proportion"` for binomial data, where they compare the predicted
#'   probability with the observed fraction of successes. The Brier score
#'   has unit `"trial"`.
#' * `training` holds the analysis `row_ids`, `coefficients` (link scale),
#'   the random-effect `covariance`, the native `dispersion` (1 for Poisson
#'   and binomial), and any fitted `scale_est`, `normalization`, and
#'   `bin_reference`.
#' * `selection_records` are always kept, with `logLik`, `K`, `native_K`,
#'   `ecological_K`, `nobs`, and the `analysis` row IDs.
#' * `diagnostics` holds `convergence`, `positive_hessian`, `native_message`,
#'   and `scale_convergence`, and `integrity_sha256` protects the records as
#'   for field results.
#'
#' # Bayesian results (multiScaleR_bayes_cv)
#'
#' [cross_validate_bayes()] scores held-out blocks with the posterior
#' predictive density, the only score in the package that averages over
#' parameter uncertainty. Its `summary` reports the joint negative log score
#' per observation, and its `predictions` hold posterior mean responses and
#' marginal row scores. It has no training likelihood, so [select_cv()] does
#' not apply. See [cross_validate_bayes()] for its elements.
#'
#' @name multiScaleR_cv-class
#' @aliases multiScaleR_choice_cv-class multiScaleR_field_cv-class
#'   multiScaleR_mixed_cv-class multiScaleR_bayes_cv-class
#' @examples
#' data("landscape_counts")
#' data("surv_pts")
#' # Three spatial folds. The simulated coordinates are in meters.
#' folds <- spatial_cv_folds(crds(vect(surv_pts)), block_size = 1000, k = 3,
#'                           crs = 32617)
#' cv <- cross_validate_glm(counts ~ site, landscape_counts, folds,
#'                          family = "poisson", keep_selection = TRUE)
#' class(cv)
#' cv$complete           # TRUE: every fold was fitted and scored
#' cv$summary            # pooled scores, one row per metric, lower is better
#' cv$metrics[cv$metrics$metric == "rmse", ]   # the same score in each fold
#' head(cv$predictions)  # one row per held-out observation, with its losses
#' cv$fold_results       # rows fitted and held out, and the status of each fold
#' cv$training$fold1$coefficients   # log-link coefficients fitted in fold 1
#'
#' # What select_cv() reads from each fold
#' cv$selection_records$fold1$data[c("logLik", "K", "nobs")]
#'
#' # The pooled RMSE is the square root of the mean squared error over all
#' # held-out rows, not the mean of the three fold RMSEs
#' sqrt(mean(cv$predictions$rmse))
#' @seealso [cross_validate()], [cv_tab()], [select_cv()], and
#'   `vignette("v08-cross-validation", package = "multiScaleR")`.
NULL
