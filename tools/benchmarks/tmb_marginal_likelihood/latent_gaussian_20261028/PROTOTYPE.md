# Private spatial-scale workflow

The new code is an internal prototype, not an exported fitting API. Existing
multiScale_bayes() formulas, defaults, classes, diagnostics and saved fits retain
their behavior. Do not coerce a prototype to the public Bayesian class to bypass
its prediction or scalar-integration restrictions.

## Representation and meaning

One landscape predictor receives a proper finite Gaussian basis field on log
kernel scale. Projected coordinate columns, knot locations and the smoothing
range are explicitly supplied. The range stays fixed in this first model.
The coefficient covariance is proper and fixed, with whitened standard-normal
latent coefficients. The scale field shares the same conditional response and
Gaussian latent densities between Bayesian and likelihood targets.

The first prototype uses an intercept-only fixed log-scale mean and optional
existing random intercepts such as year. This keeps the field a function of
location alone. Ordinary fixed scale gradients and trends remain available in
the existing public formula interface. Random slope syntax, estimated spatial
range/correlation, AR(1), intrinsic fields and new response families are not
implemented by this batch.

The basis is centered against the training intercept and normalized so its mean
pointwise prior variance on training rows is1. After multiplying by field SD,
that SD is the root mean pointwise prior SD in log-scale units. It is not the
SD of a single realized field surface, nor the ecological kernel distance.
For Bayesian fitting its half-normal prior scale is required explicitly.
Existing annual and response-effect priors are retained.

Coordinates must be projected and their units known. Plain numeric columns do
not carry a CRS, so the prototype cannot infer or verify projection metadata.
The smoothing range uses those coordinate units; kernel scales use the input's
map units. Use consistent distance units or explicitly document conversions.
All training centering/projection/normalization is frozen for new locations.

## Internal entry points

- .msr_gaussian_block():proper fixed covariance, projection and normalized
  centered density reference. Shared components support its whitened effects.
- .msr_spatial_scale_prepare():append one spatial block to a chosen kernel.
  target="posterior" requires field_sd_prior. target="marginal_likelihood"
  accepts no population prior and retains the normalized hierarchical target.
- .msr_spatial_scale_laplace():joint posterior approximation with the original
  multi-start/stationarity/curvature/competing-mode checks. This does not establish
  Gaussian approximation accuracy or interval calibration.
- .msr_spatial_scale_objective() and .msr_spatial_scale_adfun():normalized target
  and joint/latent-Laplace adapters. The likelihood adapter is not a qualified
  optimizer, limiting-model witness or population confidence-interval interface.
- .msr_spatial_scale_predict():given a prepared object, pass parameters explicitly
  for a conditional vector. Given a fitted object, its default uses joint draws
  and returns draws-by-observations. Request scale="log" for finite log-distance
  summaries. New-group policies retain coherent shared annual effects.

The private prototype class prevents accidental public summary dispatch.
Public scale prediction and scalar/one-SD fit routes reject Gaussian projection
blocks explicitly. A spatial basis generally couples sites; it cannot be sent
through the independent-year scalar factorization.

## Interpreting maps and uncertainty

Joint draws give a scale distribution at each prediction location. Use medians,
quantiles and predefined regional contrasts. Natural-scale means can be unstable
under uncertain SDs; log summaries and distance medians/quantiles are preferable.
The finite basis implies a joint covariance at new locations through its shared
coefficients. It adds no independent new-location field innovation. Extrapolation
follows the chosen basis, not a fully resolved spatial GP.

When the landscape response coefficient is zero, its scale field is absent from
the likelihood. Priors can produce smooth maps even when the data contain little
scale information. A positive equal-tail interval for field SD is not evidence
against a constant-scale null. Data support, response-only spatial confounding,
spatial holdouts and false-positive/calibration studies remain necessary.

The prepared landscape neighborhoods still bound available information.
Large inferred scales increasingly average those neighborhoods and do not recover
landscape cells that were never extracted. The two fixed pilots in
PILOT_PROTOCOL.md are illustrations of computation and recovery, not a completed
scientific validation or a general claim of nonstationarity detection.
