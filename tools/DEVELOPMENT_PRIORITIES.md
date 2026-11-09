# Development priorities

The November 8 priority update and nonlinear evaluation place supported GDR
evaluation/refitting and retained training recipes ahead of the remaining
spatial extensions. The evaluation adds two `ns4` full references and extends
the training-contract fixtures to both three- and four-coefficient bases.

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
unit contract for native refits is added in 0.7.25 below. Learned GDR component
recipes and a frozen-reference fitting alternative still need integration.

Native exponential response fields have opt-in projected coordinate descriptors
and serial raster prediction with conservative memory planning. Conformance
covers the declared Gaussian/log-Poisson subset with one shared field and
selected global ecological scales. Native mean SEs omit ecological-scale
uncertainty and observation noise. Zero-field references are not marginal
predictions, and raster batches do not supply cross-row covariance.

## Private GDR prerequisites in 0.7.24

Version 0.7.24 adds a normalized Gaussian,
finite-Poisson-cutoff GDR joint program. It passed six fixed native-reference
likelihood/score, curvature, full covariance, joint-delta density-SE and fresh
latent-start checks, plus a normal vignette-enabled package check and the
complete exact-archive installed suite. This supplies derivative prerequisites;
it does not enable `make_gdr_refitter()` for TMB, qualify new outer fits or
connect an analytic GDR scale search. Existing admission guards remain.

Version 0.7.25 separates immutable private observation preparation from
refreshed abundance designs with exact row and column identity. It retains the
original qualified seven-moment primitive. A single-callback experiment passed
its mathematical proof but remains outside the package after an interrupted
AD-memory/cold check. Numeric fingerprints are versioned and remain distinct
from training recipes, provider/DLL provenance and inference eligibility.
Final-source qualification requires all eight references, including both ns4
examples, with fresh R process boundaries and every original numerical gate.
Public TMB fitting remains guarded.

The new evaluation supplies eight independently qualified analysis fits with
72 native/numerical checks, 24 fold fits and 20 training-contract fixtures.
Those are acceptance examples for the package adapter. They are distinct from
the separate 54-record semantic ledger and from the six 0.7.24 package
fixed-point checks. The 0.7.25 context qualification adds both ns4 references.
See `GDR_JOINT_TARGET.md` for the private likelihood contract.
The local evaluation and its acceptance ledgers are retained separately from
package source.

## Explicit native-refit units in 0.7.25

The opt-in raw transformation applies to optimized kernel/metric outputs and
retains named identity transforms through prediction, effect curves, raster
projection and spatial GLM/choice/NB2/Gamma CV initial fits. Candidate
standardization stays the default. Fixed metrics keep their existing training
transforms; upstream weighted-SQ reference normalization is unchanged. Raw mode
requires native refitting, preserving analytic derivative and Bayesian target
contracts. It adds no penalty parameter or uncertainty source. This supplies
the raw-unit fitting/prediction foundation, not GDR TMB or paired-CV admission.

## Next integration steps

1. Qualify fast GDR TMB evaluation/refitting for one primary survey, Poisson
   abundance, half-normal distance detection, the established removal law and
   supported normalized Gaussian observer/abundance effects. Reuse exact private
   thinning prerequisites. Preserve native constants, independent/cold/latent
   starts, scores, curvature, full covariance and joint-delta uncertainty gates.
   Use all eight actual linear/quadratic/ns3/ns4 references for a chosen new
   implementation. Add a deliberate scalar-result contract to the optimizer;
   its existing native S4 sample checks do not admit an arbitrary scalar/list.
   Point-only evaluations must not acquire inference eligibility. Verify which
   covariance work is omitted, retaining latent curvature and full final
   native inference. Keep evaluated controls and warm-start recipe identity.
2. Extend the implemented native raw-unit contract to retained GDR component
   recipes and a separately qualified fixed-reference fitting option. Preserve
   legacy candidate standardization, Bayesian frozen references and fixed
   observation predictors while abundance shapes or ecological scales change.
   Numeric design contexts do not reconstruct learned prediction recipes.
3. Integrate retained component recipes with paired GDR held-out scores and
   reusable fold jobs. Learn scales, transformations and starts from training
   data only. Declare integration over unknown states, retain failed folds and
   verify that unused Gaussian levels integrate to one. Port all 20 supplied
   ns3/ns4 fixtures. Keep conditional means, new-group prediction and
   uncertainty sources distinct. Preserve actual training-selected shapes and
   each paired fold outcome; exploratory fold reuse is not new confirmation.
4. Add spline rank, conditioning, support and extrapolation diagnostics and
   clearly labeled response-curve scenarios. A near-zero knot is not a
   biological threshold. Natural tails are linear on the linear predictor
   scale, and fixed-complexity splines need not be monotone or single-peaked.
   Retain requested and actual native knots. The zero-heavy vireo ns4 example
   needs that diagnostic, not a silently changed knot or biological threshold.
5. Retain authenticated provenance for new spatial preparations with an explicit
   legacy/unverified path. Bind raw coordinate predictors to geometry without
   relearning transformations from prediction batches.
6. Qualify broader native covariance/family combinations, normalized correlated
   joint fitting and spatial mixed-model CV separately. Spatially varying
   ecological kernels and penalized GDR smooths need their own reference,
   identifiability, integration and uncertainty validation.

Existing priors, fixed reference normalization, diagnostics, parameter penalties
and admission guards remain. Numeric conformance does not resolve observation-
law misspecification or establish ecological recovery. Use focused regressions,
generated help, a fresh vignette-enabled check and complete installed tests for
each coherent batch. Keep checkpoints compact and preserve failed attempts.
