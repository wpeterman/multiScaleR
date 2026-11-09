# Private Gaussian GDR joint-target contract

The internal program combines native finite-cutoff Poisson GDR margin
probabilities with every normalized Gaussian intercept density. It is a
prepared-numeric derivative prerequisite, not a fitted-model, prediction,
model-selection or uncertainty adapter. `make_gdr_refitter()` retains its
original TMB rejection. Offset, ZIP and all other existing guards remain.

The initial native layout has one primary survey, complete distance/removal
margins, half-normal distance detection, a 75 m complete point circle with
distance breaks c(0,25,75), and density per square kilometer. It retains original
native matrices and their column labels. Formula construction stays with
unmarked and existing learned-term helpers, including fixed raw-cover spline
recipes and factor interactions. No second formula parser or competing basis
system is introduced. Removal fields, annual latent availability, other
families, distance keys and partial geometries need separate qualification.
Fixed Year regression columns may remain in abundance; they are not annual
latent effects.

Each explicit numeric layout includes ordered survey/removal-interval identities,
both count margins, all three design matrices, Gaussian row-to-level maps,
original circle area, annulus areas/weights, interval lengths and cutoff. Invalid
rows, missing design columns, unequal counts, unsupported blocks and altered
geometry stop. A present field retains its complete latent vector and one
log-SD, including unused levels. An absent field has neither vector nor log-SD.
Setting modes to zero or using a small SD does not create that reduced model.

For an observed total n and cutoff K, the program adds
log P(Poisson(lambda*(1-q))<=K-n) to the native infinite-sum margin identity.
K remains computational truncation, not a renormalized abundance distribution
or a penalty parameter. It retains one observed total and the native margin
constant. A paired-cell likelihood has a different observed constant and cannot
be substituted. All Gaussian densities remain normalized.

Stable half-normal annulus moments use explicit adjoints through fourth order,
including the uniform-distance limit. The original qualified seven-moment primitive and its beyond-chain derivative
guard remain unchanged. A separate single-callback experiment is not integrated:
its mathematical derivative proof does not establish conservative AD memory.
The serial R callbacks are ephemeral and
cannot be checkpointed as portable objectives or evaluated through OpenMP.
RTMB2.0 is the qualified callback layout. Bind each callback function by value;
lazy loop arguments can resolve a later binding and change every derivative.
The [RTMB adjoint documentation](https://cran.r-project.org/web/packages/RTMB/refman/RTMB.html#ADjoint)
describes the derivative contract and thread restriction. This R callback route
has overhead, so no speed gain is claimed from operation counts alone.

The private observation context validates the invariant response layout once
and constructs its constants, geometry and serial callbacks once. Its prepare
method accepts a refreshed abundance design only with the exact original
survey identities, dimensions, column names and row labels. Distance/removal
predictors and Gaussian densities are still evaluated at the current fitted
parameters. They are not frozen probability vectors. Existing one-program
calls retain their formal interface through a convenience wrapper.

Response fingerprints bind counts, observation designs, rows, geometry, K,
parameter layout, abundance layout and the relevant R implementation. Program
fingerprints additionally bind the refreshed abundance VALUES; their private
format version is 2. Rebuild the context when invariant inputs change. Prepare
a separate program and new AD tape when abundance values change. Never mutate
env$data and assume recorded constants changed. Caller edits or later prepare
calls cannot change a previous program. A numeric design is not a retained
training/basis recipe or permission to relearn prediction transformations.
Neither fingerprint supplies complete dependency/DLL provenance, a portable
fitted model or inference eligibility.

Metadata reports fixed and latent coordinate counts separately. Fixed counts
include regression and present-field log-SDs, excluding latent modes and fixed
K. They do not authorize an AIC penalty or add an unestimated ecological scale.
No coefficient SE, kernel-scale interval, new-state variation or observation
noise is supplied. A flat distance Hessian near the uniform limit remains a
failure of inference, even when the stable likelihood value is finite.

Native near-uniform arithmetic failures remain in the development evidence.
The unmodified native half-normal path exceeded its probability gate at one
large-scale reference and returned NaN at another. A stable limiting value does
not waive those failures or silently switch the user's model to uniform.
Complete-detection/zero-abundance derivative endpoints, overflow, new outer
starts, final native inference, raw-reference estimation and prediction, and
paired held-out integration remain separate gates. Bayesian priors and frozen
normalization, existing diagnostics, penalties and uncertainty distinctions
are unchanged.
