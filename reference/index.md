# Package index

## Concepts and objects

Start here. Background on kernels, fitting methods, and predictions, a
glossary of terms, and descriptions of the objects the main functions
return.

- [`multiScaleR`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-package.md)
  [`multiScaleR-package`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-package.md)
  : multiScaleR: estimate scales of effect for spatial covariates
- [`multiScaleR-glossary`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-glossary.md)
  : Glossary of terms used in multiScaleR
- [`multiScaleR-kernels`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md)
  : Kernels and scales of effect
- [`multiScaleR-fitting-methods`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md)
  : Fitting methods for estimating scales of effect
- [`multiScaleR-prediction`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md)
  : Prediction quantities and their uncertainty
- [`multiScaleR-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-class.md)
  : Fitted multiScaleR objects
- [`multiScaleR_data-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_data-class.md)
  : Prepared inputs for scale estimation
- [`multiScaleR_bayes-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md)
  : Bayesian landscape-scale fits
- [`multiScaleR_cv-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
  [`multiScaleR_choice_cv-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
  [`multiScaleR_field_cv-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
  [`multiScaleR_mixed_cv-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
  [`multiScaleR_bayes_cv-class`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
  : Cross-validation results

## Prepare inputs

Build the kernel inputs and declare how each raster covariate is
summarized before optimization.

- [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  : Prepare Kernel Inputs for Scale Optimization
- [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md)
  : Prepare Kernel Inputs from Group-Specific Rasters
- [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`print(`*`<multiScaleR_vars>`*`)`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  : Define Model Covariates from Raster Layers
- [`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md)
  : Prepare spatial covariates for native unmarked prediction
- [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md)
  [`print(`*`<msr_unmarked_audit>`*`)`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md)
  : Audit Survey Identities and Paired Distance-Removal Counts
- [`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
  [`print(`*`<multiScaleR_ram>`*`)`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
  : Estimate Memory Requirements for Multiscale Optimization

## Optimize scales of effect

Estimate the scale of effect for each covariate by maximum likelihood
and assess the fit.

- [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  : Estimate scales of effect by maximum likelihood
- [`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)
  : Profile Model Fit Across Sigma
- [`scale_inference()`](https://wpeterman.github.io/multiScaleR/reference/scale_inference.md)
  : Profile and Bootstrap Inference for an NB2 Scale
- [`diagnostics()`](https://wpeterman.github.io/multiScaleR/reference/diagnostics.md)
  : Retrieve diagnostics from multiScaleR objects
- [`vcov(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/vcov.multiScaleR.md)
  : Covariance matrices for an analytic multiscale fit

## Bayesian scale estimation

Estimate scales of effect with full posterior uncertainty.

- [`summary(`*`<multiScaleR_bayes>`*`)`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  [`print(`*`<multiScaleR_bayes>`*`)`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  : Estimate Scales of Effect with Bayesian Inference
- [`check_bayes_setup()`](https://wpeterman.github.io/multiScaleR/reference/check_bayes_setup.md)
  [`print(`*`<multiScaleR_bayes_setup>`*`)`](https://wpeterman.github.io/multiScaleR/reference/check_bayes_setup.md)
  : Check Optional Software for Bayesian Scale Models
- [`simulate_bayes_prior()`](https://wpeterman.github.io/multiScaleR/reference/simulate_bayes_prior.md)
  : Simulate Prior Predictive Responses and Landscape Scales
- [`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md)
  : Summarize Posterior Landscape Scales
- [`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md)
  : Check Sampling Precision of Fitted Bayesian Scales
- [`diagnose_bayes_probabilities()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_probabilities.md)
  : Check Monte Carlo Precision of Scale Probabilities
- [`diagnose_bayes_laplace()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_laplace.md)
  : Assess Fitted Scale Probabilities from an Uncorrected Laplace Fit
- [`plot(`*`<multiScaleR_bayes>`*`)`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR_bayes.md)
  : Plot Bayesian Landscape Scales and Group Variation
- [`simulate(`*`<multiScaleR_bayes>`*`)`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md)
  : Simulate Posterior Predictive Responses

## Summarize and compare models

- [`summary(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md)
  : Summarize a multiScaleR fit
- [`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)
  : Convert Sigma into an Effective Distance
- [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)
  : Rank models by AIC or AICc
- [`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md)
  : Rank models by BIC
- [`check_bayes_setup()`](https://wpeterman.github.io/multiScaleR/reference/check_bayes_setup.md)
  [`print(`*`<multiScaleR_bayes_setup>`*`)`](https://wpeterman.github.io/multiScaleR/reference/check_bayes_setup.md)
  : Check Optional Software for Bayesian Scale Models
- [`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
  [`print(`*`<multiScaleR_cv>`*`)`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
  : Cross-Validate Estimated Scales with Explicit Folds
- [`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md)
  [`print(`*`<multiScaleR_choice_cv>`*`)`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md)
  : Cross-Validate Conditional Choices in Complete Matched Sets
- [`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
  [`print(`*`<multiScaleR_ram>`*`)`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
  : Estimate Memory Requirements for Multiscale Optimization
- [`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  [`print(`*`<multiScaleR_folds>`*`)`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  : Create Grouped Cross-Validation Folds
- [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  [`print(`*`<multiScaleR_vars>`*`)`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  : Define Model Covariates from Raster Layers
- [`summary(`*`<multiScaleR_bayes>`*`)`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  [`print(`*`<multiScaleR_bayes>`*`)`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  : Estimate Scales of Effect with Bayesian Inference
- [`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md)
  [`print(`*`<msr_gdr_histories>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md)
  : Convert Linked Bird Records into Paired Distance-Removal Counts
- [`prepare_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_intervals.md)
  [`print(`*`<msr_gdr_intervals>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_intervals.md)
  : Add Interval-Level Covariates to a Distance-Removal Frame
- [`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
  [`predict(`*`<multiScaleR_prediction>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
  [`print(`*`<multiScaleR_prediction>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
  : Prepare a fixed set of predictions for repeated use
- [`print(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/print.multiScaleR.md)
  : Print a multiScaleR fit
- [`print(`*`<multiScaleR_data>`*`)`](https://wpeterman.github.io/multiScaleR/reference/print.multiScaleR_data.md)
  : Print prepared multiScaleR inputs
- [`print(`*`<summary_multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/print.summary_multiScaleR.md)
  : Print a summary of a multiScaleR fit
- [`sim_dat_gdr()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md)
  [`print(`*`<msr_gdr_simulation>`*`)`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md)
  : Simulate Paired Distance-Removal Surveys with Retained Truth
- [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md)
  [`print(`*`<msr_unmarked_audit>`*`)`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md)
  : Audit Survey Identities and Paired Distance-Removal Counts

## Cross-validation

- [`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
  [`print(`*`<multiScaleR_cv>`*`)`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
  : Cross-Validate Estimated Scales with Explicit Folds
- [`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md)
  : Cross-Validate an Ordinary GLM Reference Model
- [`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md)
  [`print(`*`<multiScaleR_choice_cv>`*`)`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md)
  : Cross-Validate Conditional Choices in Complete Matched Sets
- [`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md)
  : Cross-Validate a Spatial Response Field
- [`cross_validate_mixed()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_mixed.md)
  : Cross-Validate Mixed-Effects Models
- [`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md)
  : Cross-Validate Bayesian Scale Models
- [`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
  : Compare Cross-Validated Prediction Scores
- [`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
  : Validate Training-Only Model Selection Within CV Folds
- [`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  [`print(`*`<multiScaleR_folds>`*`)`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  : Create Grouped Cross-Validation Folds
- [`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md)
  : Create Spatially Blocked Cross-Validation Folds

## Prediction and spatial projection

Apply fitted scales to new observations and to rasters.

- [`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
  : Apply Kernels to Raster Layers for Mapping
- [`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
  [`predict(`*`<multiScaleR_prediction>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
  [`print(`*`<multiScaleR_prediction>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
  : Prepare a fixed set of predictions for repeated use
- [`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md)
  : Predict Landscape Scales for New Covariates and Groups
- [`predict(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)
  : Predict from a fitted multiscale model
- [`predict(`*`<multiScaleR_bayes>`*`)`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR_bayes.md)
  : Predict Means from a Bayesian Scale Model
- [`prediction_contract`](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md)
  : Read prediction contracts

## Distance-removal models

Data preparation, fitting helpers, and audits for generalized
distance-removal (GDR) models, including spatial response fields.

- [`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md)
  [`print(`*`<msr_gdr_histories>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md)
  : Convert Linked Bird Records into Paired Distance-Removal Counts
- [`prepare_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_intervals.md)
  [`print(`*`<msr_gdr_intervals>`*`)`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_intervals.md)
  : Add Interval-Level Covariates to a Distance-Removal Frame
- [`expand_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/expand_gdr_intervals.md)
  : Expand Site and Survey Covariates into GDR Intervals
- [`gdr_start_values()`](https://wpeterman.github.io/multiScaleR/reference/gdr_start_values.md)
  : List the Starting-Value Order of a Distance-Removal Model
- [`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md)
  : Create an Audited Refit Function for Distance-Removal Models
- [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md)
  : Fit a Fixed-Range Spatial Abundance Field to a Distance-Removal
  Model
- [`response_field_geometry()`](https://wpeterman.github.io/multiScaleR/reference/response_field_geometry.md)
  : Bind Projected Coordinates to a Fitted Response Field
- [`confint(`*`<unmarkedFitGDRField>`*`)`](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
  [`update(`*`<unmarkedFitGDRField>`*`)`](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
  [`simulate(`*`<unmarkedFitGDRField>`*`)`](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
  [`profile(`*`<unmarkedFitGDRField>`*`)`](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
  [`ranef(`*`<unmarkedFitGDRField>`*`)`](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
  : Fitted Distance-Removal Model with a Fixed-Range Abundance Field
- [`predict(`*`<unmarkedFitGDRField>`*`)`](https://wpeterman.github.io/multiScaleR/reference/predict-unmarkedFitGDRField-method.md)
  : Predict from a Distance-Removal Model with a Fixed-Range Abundance
  Field
- [`vcov(`*`<unmarkedFitGDRField>`*`)`](https://wpeterman.github.io/multiScaleR/reference/vcov-unmarkedFitGDRField-method.md)
  : Covariance Matrix for a Distance-Removal Model with a Fixed-Range
  Field

## Plot

- [`plot(`*`<multiScaleR>`*`)`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR.md)
  : Plot method for multiScaleR objects
- [`plot(`*`<sigma_profile>`*`)`](https://wpeterman.github.io/multiScaleR/reference/plot.sigma_profile.md)
  : Plot a Sigma Profile
- [`plot_kernel()`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md)
  : Plot a Kernel Weight Curve
- [`plot_marginal_effects()`](https://wpeterman.github.io/multiScaleR/reference/plot_marginal_effects.md)
  : Plot marginal effects from a fitted model

## Simulation

- [`sim_rast()`](https://wpeterman.github.io/multiScaleR/reference/sim_rast.md)
  : Simulate Spatially Autocorrelated Raster Surfaces

- [`sim_dat()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat.md)
  : Simulate Responses at Known Scales of Effect

- [`sim_dat_unmarked()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_unmarked.md)
  :

  Simulate data for optimizing scales of effect with `unmarked`

- [`sim_dat_gdr()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md)
  [`print(`*`<msr_gdr_simulation>`*`)`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md)
  : Simulate Paired Distance-Removal Surveys with Retained Truth

## Example data

- [`count_data`](https://wpeterman.github.io/multiScaleR/reference/count_data.md)
  : Simulated counts for the hab example
- [`landscape_counts`](https://wpeterman.github.io/multiScaleR/reference/landscape_counts.md)
  : Simulated counts for the vignette example
- [`hab`](https://wpeterman.github.io/multiScaleR/reference/hab.md) :
  Example habitat raster
- [`landscape`](https://wpeterman.github.io/multiScaleR/reference/landscape.md)
  : Example landscape raster
- [`pts`](https://wpeterman.github.io/multiScaleR/reference/pts.md) :
  Survey points for the hab example
- [`surv_pts`](https://wpeterman.github.io/multiScaleR/reference/surv_pts.md)
  : Survey points for the vignette example
