# Spatial response contracts in 0.7.30

The coordinate-aware native glmmTMB adapter now covers one shared exponential,
Gaussian or Matern field with Gaussian, log-Poisson, logit-binomial, log-NB2 and
log-Gamma responses. It retains the qualified glmmTMB 1.1.14 native likelihood,
predictions and mean-prediction SEs. Projected two-dimensional CRS, exact row
order, known units, geometry/source binding and memory/serial guards remain.
Other random effects in the same coordinate-aware map are still guarded.

Fifteen prescribed native constructions verify point, snapshot and reordered/
repeated/new coordinates. Their ranges and Matern shapes are fixed through
native parameter maps. Native logLik degrees of freedom count only estimated
parameters, including ordinary coefficients, field SD and applicable dispersion.
Realized field states and fixed range/shape add no separate penalty. The
structural covariance-coordinate diagnostic continues to include mapped
coordinates; changing that saved descriptor would invalidate old bindings.
The original fully estimated exponential references remain unchanged.

Three fixed-reference scale refits and an actual NB2/Matern raster verify the
expanded routes. The refit fallback resolves a native glmmTMB constructor in
its retained formula environment, preserving locally bound family/settings.
Ordinary update remains the first route. Its original failure is retained as
fallback metadata. Replacement data, controls, maps, priors and native counts
remain explicit. Mixed-model effect curves now attach their native zero-effect
response prediction contract without changing curve values or confidence limits.

Binomial output is the native per-trial probability, including grouped counts.
Poisson/NB2 output is an expected count, and Gaussian/Gamma output is a response
mean. Native mean SEs are conditional on selected ecological scales and add no
observation noise or ecological-scale uncertainty. Setting fields to zero is
not marginalization over them. Full raster cross-row covariance is not returned.

A separate private GDR context supplies a fixed-range exponential abundance
field to the original finite-K Poisson/half-normal distance/removal likelihood.
It uses correlation-whitened states with the original normalized Gaussian SD,
one consistent observed location per state and no annual field. Frozen range
adds zero estimated parameters. Six declared coefficient/SD points compare
joint derivatives, Laplace targets/scores, original/cold conditional modes and
positive current-mode latent precision against native compiled GDR with the
same transformed random design. These are not optimized fixed-parameter fits.

Its compact portable numeric projection supplies conditional new-field residual
variance, multiplied by exp(2*lsigma_lambda). It omits parameter uncertainty,
training-field uncertainty, ecological-scale uncertainty, observation noise
and cross-row field covariance. It is not a fitted response or mean-prediction
SE. The full likelihood context remains serial and nonportable. Public spatial
GDR fitting, prediction, AIC and CV admission remain false.

No prior, normalized target, uncertainty gate, parameter penalty or existing
independence/zero/limiting/curvature check is relaxed. These are software
compatibility and derivative qualifications, not recovery or calibration claims.
