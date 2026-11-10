# Response spatial integration priorities

Priority is spatial dependence in the response, distinct from a spatially
varying ecological kernel. The latter remains private and retains its failed
reference-precision/tail gates. A response field must retain the model's actual
observation likelihood, normalization and covariance interpretation.

## Current verified boundary

Native glmmTMB refits preserve the fitted model and native covariance. Field
metadata identifies exponential, Gaussian and Matern native structures, but
coordinate-aware raster qualification currently admits only one shared
exponential field, Gaussian or log-Poisson response, and no other random fields.
The joint analytic mixed backend remains limited to a single iid intercept.
Native zero-field predictions are not random-field integration. Mean SEs
condition on ecological scales and omit observation noise.

Native GDR TMB refits and retained component recipes are implemented. Their iid
abundance blocks are not continuous spatial fields. The guarded private GDR
joint program preserves normalized Gaussian states, the finite Poisson cutoff,
distance/removal margins, factor/spline designs and current-mode curvature.

## Next bounded batches

1. Qualify more native glmmTMB field/family combinations using genuine native
   fitting and predictions at known/new/repeated/permuted coordinates. Include
   binomial, NB2 and Gamma where appropriate, as well as Gaussian/Matern
   covariance. Retain provider-version, exact CRS/row, memory, covariance and
   explicit-new-level guards. Test serial points/raster chunks, snapshots,
   effects and native parameter counts. Failed native fits remain failures.
   Do not claim recovery or comparative speed from software fixtures.
2. Add a separate private GDR fixed-range response-field prerequisite. Bind
   projected coordinates to exact abundance-state indices and use a frozen
   correlation Cholesky to whiten the normalized Gaussian field. Keep the GDR
   observation target and all other parameter blocks unchanged. No year term
   is introduced. Verify likelihood/score/Hessian against a native TMB target
   with the corresponding transformed random design and independent density
   normalization. Check the independence limit and proper range/SD counts.
3. Develop bounded new-coordinate projection with explicit residual field
   variance, distinct from coefficient/native mean SE and ecological-scale
   uncertainty. No grid covariance or automatic variance combination. Keep
   public GDR spatial fitting/prediction guarded until optimized original/cold,
   curvature, full covariance and reference gates have passed.

## Model-class boundary

GLM/lm have no native random field; adding one changes their statistical model.
mgcv spatial smooths, spaMM spatial covariance, nlme residual correlation,
glmmTMB correlated random fields and custom integrated/Bayesian latent fields
need separate likelihood and covariance qualification. lme4 iid grouping alone
does not become continuous spatial dependence. Conditional-logistic/step-
selection spatial effects also need their matched risk-set likelihood retained.
Do not infer cross-class support from accepting a spatial formula symbol.

Every coherent source batch needs focused tests, generated docs, a fresh normal
vignette-enabled check and complete exact-archive tests. Preserve guards,
normalized targets, priors, parameter counts, failures and owner edits. No
undeclared recovery, optimizer, sampler or empirical full-N panel is authorized
by this plan. Declare any new bounded software constructions before launch.
