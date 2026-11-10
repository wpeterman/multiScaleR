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
unit contract for native refits is added in 0.7.25 below. Trained GDR prediction
recipes have their separate 0.7.28 qualification below. A frozen-reference
fitting alternative still needs integration.

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

## Native TMB refitting in 0.7.27

The native GDR refitter accepts explicit evaluated TMB settings and fixed/log-SD
constructor starts, retains all original/cold attempts, and returns the actual
native S4 model. The native compiled likelihood and Gaussian densities remain
unchanged. Controls are forwarded locally, and native sdreport is evaluated at
the optimizer coordinates. Intermediate TMB attempts retain full covariance
checks; no likelihood-only efficiency claim follows. TMB ZIP, offsets and
unusable native covariance remain failures. Priors, counts, fixed references
and kernel-conditional uncertainty conventions are unchanged.

## Native GDR prediction recipes in 0.7.28

Component-scoped C/TMB GDR predictions retain checked original fixed designs,
factor coding and supported learned transformations. Site, primary-period and
removal-interval recipes remain distinct through points, frames, numerical
snapshots, raster chunks and active-component curves. Known-group columns and
native covariance are retained; zero effects are not random-effect integration.
Original covariates are reused within chunk contexts. Conditional learned
random designs, offsets, TMB ZIP, missing coding and changed layouts remain
guarded. No new prior, target, penalty parameter, paired-CV admission or scale-
uncertainty claim follows.

## Next integration steps

The November 9 additions put faster development feedback ahead of the next
package batch. The validation tools now inventory every test file and provide
explicit partial feedback profiles, with compact block/file timing records for
future runs. The normal vignette-enabled check and complete exact-archive suite
remain publication gates. Use measurements to remove redundant fixture work;
do not trim scientific guards or use assertion totals as a growth target.

The added OBBA_lsm profile-GLM memo concerns 14,347 binary recorded-detection
responses, a different target and cohort from paired GDR surveys. It requests
an opt-in raw/reference analytic transformation contract with fixed physical-
unit ns2 recipes, then exact native hard-radius MI/ED search outside smooth
bandwidth/coefficient profiling. Keep hard-radius differentiability guards
until that hybrid route is separately validated. Preserve exact event/edge
inclusion, finite-coefficient existence checks, failed/unknown/lower points,
boundary and multiple-basin evidence, and radius-versus-smooth uncertainty.
Nominal coefficient-plus-scale counts are 7/10 for its declared amount/metric
families; unfinished searches do not supply maximized AICc.

A prepared-feature provider must represent actual immutable inputs and site/
response identities without fabricated cell lists. Bound exact-radius caches,
qualify cached/uncached and fresh-process equality, and retain prediction
recipes. Benchmark equal scientific work only after qualification and when
whole process-tree memory permits independent species. Do not copy its large
feature store or change pinned analysis libraries. These are new priorities,
not implemented hybrid fitting, public prepared-provider admission or a speedup
claim. Retained GDR prediction recipes below remain a correctness priority.

1. Qualify a faster point-only GDR TMB interface after native refitting. Start with one primary survey, Poisson
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
2. Extend the implemented native raw-unit contract with a separately qualified
   fixed-reference fitting option, retaining the checked GDR prediction recipes. Preserve
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
4. Version 0.7.26 adds opt-in rank/conditioning summaries and exactly verified
   state spline recipes to the existing fixed-effect native C GDR audit domain.
   It retains requested calls versus actual knots, native warnings/failures and
   compact training/scenario support counts, without fitting or prediction
   admission. Broader TMB/random-effect component diagnostics and response-curve
   scenarios remain separate. A near-zero knot is not a
   biological threshold. Natural-spline basis columns have linear tails;
   added products or transformations can change the full predictor's tails.
   Fixed-complexity spline responses need not be monotone or single-peaked.
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

## Updated profile-GLM memo queue, November 9

[The ordered implementation plan](benchmarks/tmb_marginal_likelihood/profile_glm_integration_20261109/NEXT_STEPS.md) incorporates the memo's 15:26/16:06 UTC updates and supplements the priorities above. After the GDR prediction batch, qualify an official prepared raw-design binomial evaluator and independent finite-MLE status, then opt-in physical spline references, exact hard-radius hybrid search, bounded authenticated providers/caches, downstream prediction/selection/CV and equal-work efficiency. External fixed-point benchmarks are attributed evidence to port, not new package acceptance or whole-optimizer speedup claims. Preserve the current hard-radius guard, exact sample/response identity, nominal 7/10 counts, failures and all uncertainty distinctions.

## Additional GDR follow-up, November 9

[The updated GDR queue](benchmarks/tmb_marginal_likelihood/pb_20261109/GDR_NEXT16.md)
incorporates the OBBA2 native-TMB evaluation memo. The native adapter and
0.7.28 component recipes are implemented. Remaining GDR work is qualified
point-only finite-K refitting with explicit numeric result contracts, paired
prediction/fold scoring, scaled convergence and response-encoding metadata,
and a supported controlled initial-fit route. Penalized smooths and continuous
spatial abundance remain separate evidence-led extensions. Preserve all32
full/fold cases,20 training/unused-group fixtures and the analysis-driver
factorial-constant correction. These join the current raw-binomial/provider
priorities; they do not replace the GDR observation law or loosen admission.

## Prepared binomial and current spatial priority, 0.7.29

[Prepared raw-binomial contracts](benchmarks/tmb_marginal_likelihood/pb_20261109/CONTRACT16.md) now have exact original-design existence certificates, focused installed/worker tests and full package qualification. Public raw fitting/selection/provider and uncertainty guards remain closed. The latest direct priority is [response spatial integration](benchmarks/tmb_marginal_likelihood/pb_20261109/SPATIAL_PLAN.md), especially native glmmTMB and a proper normalized GDR abundance field. Existing iid groups and spatially varying ecological kernels remain distinct.
