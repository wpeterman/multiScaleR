# Evidence for experimental and approximate inference paths

The package uses several statistical and numerical targets. The initial
qualification and subsequent studies retain failed fits and unavailable
comparisons. The links below provide the settings and denominators behind
each result.

| Workflow | Evidence and decision |
|:--|:--|
| Analytic GLM fitting and the NB2 Poisson face | The [independent likelihood study](RESULTS.md) covers 120 fits. Native and analytic likelihoods agree when they use the same dispersion range. Point estimation is supported for the tested independent-response models. Wald coverage remains experimental. |
| Analytic Wald intervals | The [2,016-dataset study](../analytic_wald_intervals/REPORT.md) compares reported Wald coverage with independent likelihood profiling across families, sample sizes, and identification scenarios. Matched local curvature can agree while Wald coverage remains below nominal. |
| Joint mixed scale fitting | The [joint mixed assessment](RESULTS.md) checks native likelihoods across Poisson, NB2, and zero-inflated Poisson random-intercept simulations. The fitted scale is a point estimate. Native standard errors condition on that scale. |
| NB2 likelihood profile and bootstrap | The [pilot](RESULTS.md) checked the rejection cutoff. The [replicated pointwise study](../nb2_bootstrap_replicated/RESULTS.md) completed 1,176 paired refits. A [larger pointwise study](../nb2_bootstrap_grid_sets_20261113/POINTWISE_RESULTS.md) compared bootstrap and chi-square decisions at the known scale. The [48-dataset grid study](../nb2_bootstrap_grid_sets_20261113/GRID_RESULTS.md) found that coverage depended on the grid even when shared-node decisions were held fixed. Calibration remains experimental. |
| Hard-radius landscape profiles | The [replicated radius-profile assessment](../radius_profile_20261112/REPORT.md) reproduced the independently selected grid row in all 12 Gaussian and conditional-choice datasets. Flat and boundary fits returned point estimates with unavailable scale SEs. The assessment addresses likelihood agreement and selection on the declared grid. |
| Smooth conditional-logistic bandwidths | The [one-case choice study](../conditional_choice_bandwidth_recovery/REPORT.md) scored all 384 simulated datasets. Weak habitat effects could improve held-out choice prediction while leaving the bandwidth poorly identified. |
| Native Cox transformation metadata | The [provider compatibility check](COX-PROVIDER.md) identifies differences in retained transformation metadata and verifies conditional-logistic prediction and refits with two survival versions. |
| Optional unweighted surface roughness cache | The [paired scalar comparison](../surface_sq_cache_20261112/REPORT.md) reproduced all 24 declared-grid selections. A large elevation datum changed four continuous radius estimates by at most 2.94e-10 map units. The cache remains opt-in. |
| Mixed-effects CV | The [replicated mixed CV assessment](../mixed_cv_replicated/REPORT.md) compares exact Gaussian conditioning and independent GLMM quadrature. It supports the fitted-parameter predictions and scores for eligible models. Native fit availability varies by family and random-effect structure. |
| Spatial field CV | The [quadrature assessment](../spatial_field_cv_quadrature/REPORT.md) covers 1,632 scenarios. The worst non-Gaussian joint score error was 0.0483 nats per assessment observation. Tilted Laplace moments improved log-link predictions. A [follow-up with 12 and 22 field locations](../spatial_field_cv_larger/REPORT.md) compared 28 approximate scores and predictive moments with independent posterior simulation. The [fitted-parameter study](../spatial_field_cv_fitted/REPORT.md) scored 341 of 360 independent fields. Plug-in score and variance reliability depended on field geometry and parameter estimation. |
| Bayesian CV | The [new-group integration assessment](../bayesian_scales/cv_new_group_integration/REPORT.md) checks simulation of shared new effects against analytic Gaussian and independent Poisson references. A [refitted Gaussian and NB2 study](../bayesian_scales/cv_refits_20261113/REPORT.md) scored all 64 training refits and found three close model comparisons reversed by the default integration budget. Rare joint responses can need more new-effect draws or greater posterior precision. |
| Full-posterior slice sampling | A [matched slice and NUTS assessment](../bayesian_scales/uncorrected_laplace_20261112/REPORT.md) returned four eligible Gaussian and NB2 fits with zero material differences among 28 fitted-scale probability comparisons. Existing tests also compare slice sampling with Gaussian quadrature. |
| Uncorrected Bayesian Laplace fits | The [matched simulation](../bayesian_scales/uncorrected_laplace_20261112/REPORT.md) found material differences from NUTS in 49 of 87 joint and 39 of 71 marginal annual-scale probability comparisons. A [35-dataset follow-up](../bayesian_scales/laplace_adequacy_20261113/REPORT.md) qualified an optional discrepancy check for individual fitted scale probabilities. Most checks were unresolved. Both approximations remain exploratory. |
| Importance-corrected marginal Bayesian fitting | The [prior-predictive assessment](../bayesian_scales/qualification_20261112/RESULTS.md) contains 1,152 fits across eight designs. Both proposals met the availability criterion in seven designs. The conditional coefficient proposal also met it for repeated-site NB2. Fitted-scale probabilities, intervals and joint predictive scores have distinct evidence. |

The reference and streaming smooth-kernel implementations have direct value,
gradient and finite-difference tests in `test-smooth-kernel-engine.R` and
`test-kernel-stream.R`. Those tests support the numerical calculation.
They do not supply independent evidence about interval coverage.

## Private reference implementations

The private spatial-scale, marginal-likelihood, and factorized-likelihood
components support numerical audits and development comparisons. Their
mathematical checks include normalized density references, derivative
comparisons, zero-variance limits, and shared latent-state dependence.
Relevant tests include `test-latent-gaussian.R`,
`test-factorized-scale-likelihood.R`, `test-marginal-scale-likelihood.R`,
`test-spatial-scale-prototype.R`, `test-spatial-scale-response.R`, and
`test-spatial-response-log-lik.R`. These components have limited or private
admission in the capability catalog. Evidence for prediction conditional on
a frozen field or population estimate does not establish calibrated scale
uncertainty.

For GDR, `test-gdr-exact-primitives.R`, `test-gdr-finite-k.R`,
`test-gdr-zero-sd-parent.R`, and `test-gdr-response-field-fit.R` check
normalized reference likelihoods, limiting cases, and comparison with native
fits. The public fixed-range GDR interface has its own audited provider and
data contracts. Unknown-range integration and paired-process extensions
remain research references.

## What to establish for a broader claim

Use an independently implemented likelihood, predictive density, or known
limiting model to test the numerical target. Vary the features that make the
approximation difficult, including weak information, boundary solutions,
dependence, and rare outcomes. Record all planned fits and whether each
reference is available. Report interval coverage separately from numerical
eligibility and distinguish prediction error from Monte Carlo error.
