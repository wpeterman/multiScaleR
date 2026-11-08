# Development priorities

The November 8 priority update places GDR evaluation/refitting and retained
training recipes ahead of the remaining spatial extensions.

## Implemented in 0.7.23

Fixed-knot cubic `splines::ns()` and `splines::bs()` responses update
their native basis values and complete kernel derivatives during supported
analytic scale searches. All basis columns share the original source's
ecological scale. Factor interactions, retained prediction recipes, coefficient
SEs and parameter counts keep their existing contracts. Adaptive quantile
knots, repeated/boundary knots and other degrees require native refitting.
These are unpenalized bases, with no smoothing-variance parameter.

Knots are fixed in the modeled predictor's units. Analytic GLM, conditional-
logistic and joint mixed methods retain candidate-specific standardization
with complete derivatives. Those standardized knots do not fix raw-cover
thresholds. Bayesian targets keep the frozen reference transformation. A raw-
unit/fixed-reference GDR context remains part of the next integration work.

Native exponential response fields have opt-in projected coordinate descriptors
and serial raster prediction with conservative memory planning. Conformance
covers the declared Gaussian/log-Poisson subset with one shared field and
selected global ecological scales. Native mean SEs omit ecological-scale
uncertainty and observation noise. Zero-field references are not marginal
predictions, and raster batches do not supply cross-row covariance.

## Next integration steps

1. Qualify fast GDR TMB evaluation/refitting for one primary survey, Poisson
   abundance, half-normal distance detection, the established removal law and
   supported normalized Gaussian observer/abundance effects. Reuse exact private
   thinning prerequisites. Preserve native constants, independent/cold/latent
   starts, scores, curvature, full covariance and joint-delta uncertainty gates.
   Point-only evaluations must not acquire inference eligibility.
2. Integrate retained component recipes with paired GDR held-out scores and
   reusable fold jobs. Learn scales, transformations and starts from training
   data only. Declare integration over unknown states, retain failed folds and
   verify that unused Gaussian levels integrate to one. Keep conditional means,
   new-group prediction and uncertainty sources distinct.
3. Add spline rank, conditioning, support and extrapolation diagnostics and
   clearly labeled response-curve scenarios. A near-zero knot is not a
   biological threshold. Natural tails are linear on the linear predictor
   scale, and fixed-complexity splines need not be monotone or single-peaked.
4. Retain authenticated provenance for new spatial preparations with an explicit
   legacy/unverified path. Bind raw coordinate predictors to geometry without
   relearning transformations from prediction batches.
5. Qualify broader native covariance/family combinations, normalized correlated
   joint fitting and spatial mixed-model CV separately. Spatially varying
   ecological kernels and penalized GDR smooths need their own reference,
   identifiability, integration and uncertainty validation.

Existing priors, fixed reference normalization, diagnostics, parameter penalties
and admission guards remain. Numeric conformance does not resolve observation-
law misspecification or establish ecological recovery. Use focused regressions,
generated help, a fresh vignette-enabled check and complete installed tests for
each coherent batch. Keep checkpoints compact and preserve failed attempts.
