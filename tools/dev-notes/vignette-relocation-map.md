# Vignette relocation map

Working ledger for the documentation restructuring on branch
`docs-restructure`. It records where each block removed from a vignette goes,
so no user-relevant content is lost. Delete this file when every row reads
"done".

Line numbers refer to the file as of commit `bfdc133e3` (multiScaleR 0.7.36).
Recover the original text with, for example,
`git show bfdc133e3:vignettes/inference-contracts.Rmd`.

## `vignettes/inference-contracts.Rmd`

Rewritten as "What Your Estimates and Intervals Mean" (929 lines to about 345).

| Original lines | Content | Destination | Status |
|:---------------|:--------|:------------|:-------|
| 19-29 | Opening and prerequisites | Rewritten in place (section 1) | done |
| 31-38 | Dot formulas, removed variables, and retained transformations in `lm`, `glm`, and Cox predictions | `?predict.multiScaleR` details | pending, help-page phase |
| 39-48 | `plot_marginal_effects(type = "det" / "state" / "lambda")` for `occu()` and `pcount()` fits | `?plot_marginal_effects` | pending, help-page phase |
| 50-71 | Likelihood, marginal likelihood, and posterior targets; workflow table | Rewritten in place (section 5) | done |
| 73-77 | Hierarchical marginal-likelihood fitting is development work | `tools/dev-notes/validation-status.md` | done |
| 79-88 | Analytic backend table | Simplified table kept (section 2) | done |
| 90-113 | Weighted SQ requirements, fixed-radius metrics, `covariate_transform = "raw"`, shared preparation | `feature-compatibility` vignette and the fitting-methods help topic | vignette half done (`feature-compatibility`); fitting-methods help topic pending |
| 115-153 | Binary and grouped-binomial responses, weights, comparison with `glm()` | User Guide, binomial occurrence section, and `?multiScale_optim` | vignette half done (User Guide, "Model binary occurrence"); `?multiScale_optim` pending |
| 155-161 | Integrated likelihoods versus fixed-predictor GLMs | First four sentences folded into section 5; last sentence to `validation-status.md` | done |
| 163-191 | Parameter counting and used observations for `unmarked::pcount` comparisons; `K` as an integration bound | User Guide, `unmarked` section. One sentence kept in section 5 | partial: K and `pcount` in the User Guide; random-effect SD counting and effective-sample-size caveats pending for `?aic_tab` |
| 193-237 | Step-selection models with `backend = "conditional_logit"` and `amt` | New step-selection vignette | done in `step-selection`, except `amt` movement corrections, AICc/BIC event counts, and factor interactions (not verified, not carried) |
| 239-271 | Reading a prepared Bayesian specification (`method = "prepare"`) | `bayesian-scales` vignette | done |
| 273-317 | Expected response versus observation replicate | Rewritten in place with `simulate()` added (section 4). Full treatment stays in `bayesian-prediction` | done |
| 319-333 | Known and unseen groups, response versus scale year effects | `bayesian-prediction` vignette (already its owner) | done (offsets from this block not mentioned) |
| 335-348 | Uncertainty by prediction route; `re.form = NA`; unmarked components | Table merged into sections 2 and 6. `re.form` and unmarked component detail to `?predict.multiScaleR` | partial |
| 350-356 | `prepare_prediction()` snapshots | `?prepare_prediction` | pending, help-page phase |
| 358-388 | Units, support, numerical precision | Rewritten in place (sections 7 and 8) | done |
| 390-402 | Quick reference | Replaced with a code block, function table, and See also | done |
| 404-479 | `cross_validate()` workflow, scores, nuisance parameters | `cross-validation` vignette and `?cross_validate` | partial: workflow and scores done; Gamma dispersion and nuisance detail pending for `?cross_validate`, `?cross_validate_glm` |
| 481-517 | Grouped and spatial folds, `select_cv()`, `cross_validate_clogit()` | `cross-validation` vignette | done |
| 519-537 | Private likelihood prediction components | `validation-status.md` | done |
| 539-554 | Recovery of learned transformations in native predictions | `?predict.multiScaleR` | pending, help-page phase |
| 556-643 | glmmTMB spatial response fields, point and raster prediction | New spatial-fields vignette | done in `spatial-response-fields`, except the `response_field_geometry()` raster route (fails on real fits; see code issues) |
| 645-694 | Unmarked effect curves, raw reference rows, GAM and Cox notes | `?plot_marginal_effects`, `?prepare_prediction` | pending, help-page phase |
| 696-720 | Component-specific unmarked prediction preparation | `?prepare_prediction` | pending, help-page phase |
| 722-752 | Private GDR response-field references and correlation-range prerequisites | `validation-status.md` | done |
| 754-839 | Native fixed-range GDR abundance fields (`fit_gdr_response_field()`) | Distance-removal vignette (user part) and `?fit_gdr_response_field` | partial: user part in the distance-removal vignette; field map step fails; covariance detail pending for `?fit_gdr_response_field` |
| 841-858 | Gaussian response-field cross-validation | `cross-validation` vignette and `?cross_validate_field` | done as a short section; detail pending for `?cross_validate_field` |
| 860-900 | GDR optimizer points and `tmb_settings` | `?make_gdr_refitter` | pending, help-page phase |
| 902-929 | Remaining development gates; private finite-cutoff reference | `validation-status.md` | done |

## `vignettes/quickstart.Rmd`

Rewritten in place (395 lines to about 370; prose shorter, interpretation and
tables added).

| Original lines | Content | Destination | Status |
|:---------------|:--------|:------------|:-------|
| 57-66, 72-74 | Pointers to five other vignettes before any code | "Where to go next" table at the end | done |
| 130-133 | `kernel_prep_by_group()` pointer | "Where to go next" table | done |
| 174-192 | Parallel RAM guard, lean mode, `estimate_multiscale_ram()` | One sentence kept. Detail belongs in `?estimate_multiscale_ram` and `?multiScale_optim` (already documented there) | done |
| 203-223 | Screened start strategy | One-sentence pointer kept. Full text stays in the User Guide | done |
| 247-261, 267-272 | Planar (2D) versus line-based (1D) distance reporting and version history | Two sentences kept. Formulas in `?kernel_dist`; version history in `NEWS` | done |
| 289-296 | Parallel profiling advice | `?profile_sigma` (already documented there) | done |

## Rollout of the remaining vignettes (2026-11-11)

Development-status passages from each vignette were appended verbatim to
`validation-status.md` under a heading per source file. The blocks below left
the vignettes and still need a home on a help page. Line numbers refer to each
file at `bfdc133e3`.

| Source vignette | Original lines | Content | Target |
|:----------------|:---------------|:--------|:-------|
| `multiScaleR_Guide` | 340-342 | RAM report detail | `?estimate_multiscale_ram` |
| `multiScaleR_Guide` | 396 | Stepped-likelihood intervals | `?profile_sigma` |
| `multiScaleR_Guide` | 465-472 | `dimension = "1d"` switch | `?kernel_dist` |
| `multiScaleR_Guide` | 744-769 | Stacked site-year `join_by` layout | `?multiScale_optim` |
| `multiScaleR_Guide` | 779-781 | `unmarked` marginal-effect conditioning | `?plot_marginal_effects` |
| `multiScaleR_Guide` | 834-898 | `prepare_unmarked_data()`, `prepare_prediction()`, `include_fixed` | `?prepare_prediction`, `?prepare_unmarked_data` |
| `multiScaleR_Guide` | 926-956 | `match_buffer` detail | `?kernel_scale.raster` |
| `multiScaleR_Guide` | 1178-1251 | `sim_rast()` and `sim_dat()` walkthrough | `?sim_dat` |
| `feature-compatibility` | 31-37, 268-319, 508-546 | Dot formulas, points and scenarios, retained transformations | `?predict.multiScaleR`, `?prepare_prediction` |
| `feature-compatibility` | 39-48, 464-487, 555-566, 605-698 | Unmarked, Cox, GAM effect curves, `reference_data` | `?plot_marginal_effects` |
| `feature-compatibility` | 330-344 | K counting and comparability checks | `?aic_tab` |
| `feature-compatibility` | 160-172 | SQ definition, expow shape for hard radii | `?msr_vars`, `?multiScale_optim` |
| `feature-compatibility` | 240-266 | FFT versus focal, `match_buffer` | `?kernel_scale.raster` |
| `cross-validation` | 114-140 | Learned-transformation recovery | `?cross_validate`, `?predict.multiScaleR` |
| `cross-validation` | 166-171, 231-258 | Binomial forms, Brier weighting, Gamma dispersion | `?cross_validate_glm` |
| `cross-validation` | 277-285, 355-381, 390-416 | NB/Gamma spatial template, raw transform, grouped inputs | `?cross_validate` |
| `cross-validation` | 483-505 | clogit fold advice, dot formulas, ties | `?cross_validate_clogit` |
| `cross-validation` | 516-552 | Field CV template and score | `?cross_validate_field` |
| `cross-validation` | 577-600 | `keep_fits`, `keep_selection` records | `?select_cv` |
| `distance-removal-audits` | 27-38, 67-89 | Event versus site identity, effort, missing margins | `?validate_unmarked_data` |
| `distance-removal-audits` | 104-123 | Print behavior, simulator scope | `?sim_dat_gdr` |
| `distance-removal-audits` | 191-204 | `first_cells`, `native_cells`, `native_frame` | `?prepare_gdr_histories` |
| `distance-removal-audits` | 247-272 | Expansion details, `yearly_covariates` | `?expand_gdr_intervals` |
| `distance-removal-audits` | 430-490 | Refitter contract and TMB settings | `?make_gdr_refitter` |
| `distance-removal-audits` | 580-626 | Component training recipes | `?prepare_prediction`, `?plot_marginal_effects` |
| `bayesian-scales` | 193-228 | Formula rules, splines, dot formulas | `?multiScale_bayes` |
| `bayesian-scales` | 253-260 | Bayesian fits in effect plots | `?plot_marginal_effects` |
| `bayesian-scales` | 230-283 | Weighted SQ preparation and memory | `?msr_vars` |
| `bayesian-prediction` | 97-109, 378-409 | `row_id` matching, resampling and batching | `?predict.multiScaleR_bayes` |
| `bayesian-prediction` | 148-170, 274-299 | Support columns, coverage of the 95% radius | `?predict_scales`, `?posterior_scales` |
| `bayesian-prediction` | 343-374 | Seed reuse, grouped binomial `trials` | `?simulate.multiScaleR_bayes` |
| `bayesian-computation` | 367-398 | Weighted estimator, `importance_min` | `?posterior_scales` |
| `landscape_metric_covariates` | 80-153 | Per-metric formulas | `?msr_vars` |
| `surface_metric_covariates` | 80-87 | Metric formulas | `?msr_vars` |

Removed without a new home, for the maintainer to confirm: the biased-sampling
simulation and its inside-versus-outside RMSE comparison
(`spatial_projection_clamping` 86-228, 424-469).

## Code issues found during the rollout

Reported by the rewriting agents, then reproduced independently on
2026-11-11 against the branch build (0.7.36 plus the caption change) with
`verify_issues.R` on the running example unless noted. Status as of
2026-11-12 (0.7.37, uncommitted) is in the last column. Tests named there are
in `tests/testthat/`.

| Issue | Verdict | Evidence | Status |
|:------|:--------|:---------|:-------|
| `response_field_geometry()` rejects real glmmTMB field fits | Confirmed | `R/response_field_geometry.R:117` requires `object$backend == "refit"`. A default-backend fit has no `backend` element. Setting it by hand makes the call succeed. | Fixed. A missing `backend` element is the refit backend, in both the glmmTMB and GDR geometry routes. `test-fix-field-geometry.R` |
| Field fit from default glmmTMB starts returns the starting sigmas | Confirmed | Starting model has `pdHess = FALSE`. `multiScale_optim()` returns 415 and 415 with `NA` SEs. Only glmmTMB's own convergence warning appears. | Fixed (warning). glmmTMB reports an NA log-likelihood when its Hessian is not positive definite, so every refit returned the penalty and the gradient was zero. `multiScale_optim()` now warns, names explicit `start` values, and records `start_stall`. `test-fix-optimizer-start.R` |
| `predict()` at new field positions without `allow.new.levels = TRUE` | Confirmed | Stops with "subscript out of bounds". Works with the argument. | Fixed. Clear message naming `allow.new.levels = TRUE` and `re.form = NA`. `test-fix-field-geometry.R` |
| Fixed-buffer kernel stalls from the default start | Confirmed | Default start: 415 and 415, logLik -175.04. `start_strategy = "screen"`: 410 and 724, logLik -154.63. No warning from the default run. | Fixed (warning). Every scale is a hard radius, so the finite-difference gradient is exactly zero and L-BFGS-B stops after one evaluation. Warning plus `start_stall` diagnostic. Default search unchanged (decision for Bill). `test-fix-optimizer-start.R` |
| `profile_sigma()` defaults fail on a `profile_glm` fit | Confirmed | "profile_glm sigma values must remain within the fitted bounds." Works with `sigma_range = c(100, 600)` and on a default `refit` fit. | Fixed. `exp(log(min_D))` rounded below `min_D`. Grid endpoints are now exact. `test-fix-profile-backends.R` |
| Profiling fails untidily on `conditional_logit` fits | Confirmed | `profile_sigma()` errors about model data after a `length.out` warning. `summary(profile = TRUE)` says the object lacks information. Sigma SE is `NA`. | Fixed. `profile_sigma()` profiles with native clogit refits. `summary(profile = TRUE)` stops early with one message. `test-fix-profile-backends.R` |
| Fixed-knot splines with `backend = "conditional_logit"` | Confirmed | `ns()` and `splines::bs()` with literal knots both stop with "Native dynamic spline frame does not match its fixed basis recipe." The help lists them as supported. | Fixed. `coxph()` runs its frame through `na.omit()`, which strips the basis attributes, so the basis is checked against the saved prediction call. `test-fix-clogit-splines.R` |
| `profile_glm` rejects `lm()` | Confirmed, documentation only | Clear error asking for a `glm` fit. `glm(family = gaussian())` works. The help wording "Gaussian" should say `glm()`. | Fixed. Help says to use `glm(family = gaussian())` |
| `select_cv()` on mixed Poisson and NB2 candidates | Confirmed | Poisson default metrics include `deviance`, NB2 metrics do not. Error: "Distribution CV supports rmse, mae and neg_log_score." Works when the Poisson run sets the three shared metrics. | Fixed. `test-fix-select-cv-mixed.R` |
| `?cross_validate` wrongly requires `store_cell_data = TRUE` | Rejected | `TRUE` is the `kernel_prep()` default, and `store_cell_data = FALSE` inputs fail with a clear message. The help text is accurate. | No change needed |
| `base` has no effect on `msidi` | Confirmed | Identical values (0.7305) with the default and `base = 2`. `shdi` changes (0.838 to 1.209). Either the help or the code is wrong. | Fixed. Warning for `msidi`, `msiei`, and every other metric without a logarithm. `test-fix-landscape-base-rpr.R` |
| `rpr` without `classes_max` | Confirmed | "covariate `m` is all NA", where the help says the observed class count is used. | Fixed. Immediate error. `test-fix-landscape-base-rpr.R` |
| `plot_kernel()` line red for `expow`, black otherwise | Rejected | The dashed line is red for all four kernels. | No change needed |
| User-name paths in `bayesian-annual-example.rds` | Confirmed | A names attribute holds 64 absolute paths under a user folder. The other shipped `.rds` files have none. | Fixed. `test-fix-shipped-data-paths.R` |
| Two-class "three-class" layer in the cache generator | Confirmed | `land2` ranges 0 to 1 and `data-raw/generate_vignette_cache.R` cuts at -0.5 and 0.5. | Fixed in the script. Cache not regenerated |
| Un-exporting `gdr_start_values()` strands two help pages | Confirmed | `man/fit_gdr_response_field.Rd` and `man/make_gdr_refitter.Rd` refer to it. | Resolved. Bill kept it exported (2026-11-12) |
| GDR field route returns a sigma SE of 3.4e-30, and its map step fails | Not rerun | The reported run took about seven minutes. Still an agent report only. | Rerun 2026-11-12 (about 4 min). SE of 3.4e-30 not reproduced. The search returned its default start (415 m) after 3 evaluations with a negative Hessian (SE withheld), now flagged by `start_stall`. The map step failed on the geometry `backend` check and now works |

## Still to consolidate

Five per-vignette cache scripts (`data-raw/cache_*.R`) and their
`inst/extdata/*-cache.rds` files (about 0.7 MB together) need folding into the
single regeneration script decided earlier.
