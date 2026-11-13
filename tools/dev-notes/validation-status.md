# Validation status and development notes

This file holds development-status text that used to sit in the vignettes. It
records what has been validated, what is still restricted, and why. It is for
maintainers. Users see the corresponding restriction once, in a "Supported /
Not supported" table in the relevant vignette or help page.

The package build ignores this directory (`^tools/dev-notes$` in
`.Rbuildignore`), so nothing here ships to CRAN or the website.

Each section below quotes the original vignette text verbatim and names its
source, so the wording that guided development is not lost.

## From `vignettes/inference-contracts.Rmd` (as of commit bfdc133e3)

### Hierarchical marginal-likelihood fitting (original lines 73-77)

Unpenalized hierarchical marginal-likelihood fitting and boundary-aware SD
intervals remain development work. The public Bayesian method names above do
not expose those private prototypes as likelihood fits or confidence procedures.
See `?multiScale_bayes` for each public method's supported model structures and
diagnostic requirements.

### Integrated likelihoods versus fixed-predictor GLMs (original lines 155-161)

Integrating over random scale or response effects defines a different
likelihood from a fixed-predictor GLM. A Bayesian posterior also retains its
population priors, even if no random effects are included. Shared Bernoulli
density code does not make these objectives or their intervals interchangeable.
Compare likelihood-based information criteria only for compatible likelihood
targets, observed data and parameter-count conventions. The private unpenalized
integrated-likelihood prototypes retain their development and admission guards.

### Private likelihood prediction components (original lines 519-537)

The private marginal and factorized likelihood prototypes have conditional point
snapshot components. They evaluate an explicit population and realized-effect
state at training or new locations, then reuse the computed link, expected
response and ecological-scale outputs. New locations require the original
prepared training inputs, so their bin values, counts, distances, source
declarations and retained normalization can be checked. Snapshots keep computed
outputs and compact provenance instead of a fitted model or AD tape.

These components do not widen public fitting admission. Fitted adapters require
the original successful numerical diagnostics, selected eligible start and
zero-SD reduction. Integrated effects omitted by a factorized fit require explicit
reconstructed coordinates. An omitted effect cannot silently become zero.
Explicit zero-effect policies describe a reference prediction. The output does
not integrate population-parameter uncertainty, conditional latent variation,
new-group variation or observation noise, and supplies no SE or interval.
Numerical prediction validity does not establish uncertainty calibration.


### Private GDR response-field references (original lines 722-735)

A fixed-range GDR abundance field is an internal prerequisite, not a public
fitted-model route. It preserves the finite-cutoff Poisson/half-normal
distance/removal target, its observation constants and normalized Gaussian SD.
The canonical two- and three-distance-bin references retain a complete
75-meter sampled circle. Joint derivatives, conditional modes and curvature
are distinct from optimized fixed-parameter covariance. Optimized reference
checks use the existing 160-survey mixed fixture, both declared coefficient/SD
starts and full fixed/latent cross-covariance. Conditional log-density
contrasts retain those cross-terms; they do not establish response-mean
integration, interval calibration, spatial CV or public AIC admission. New-field
residual variance remains separate from fitted-field/parameter uncertainty,
observation noise and ecological-scale uncertainty. No year field is added.

### Private correlation-range prerequisites (original lines 737-752)

A private GDR target adds one log correlation-range parameter to the same
finite-K Poisson/half-normal observation law. Abundance states use physical
field coordinates. Their correlated Gaussian density retains its determinant
and original field SD; the range parameter adds no prior or fitting-coordinate
Jacobian. Physical and whitened joint densities differ by a state-coordinate
Jacobian, which cancels when those states are integrated. Independent stationary
AR1 and compiled native sections check these distinctions, including range/SD/
state derivatives and the existing conditional-mode curvature checks.

This is a likelihood prerequisite, not an estimated-range fitted-model route.
Its extra coordinate does not enable AIC, spatial CV, new prediction uncertainty
or public RTMB fitting. The published native GDR field still requires a supplied
fixed range. Estimated range needs independent optimum/identifiability checks,
correct penalties and full prediction covariance before its guards can change.

### Remaining development gates for native GDR fits (original lines 902-906)

Final native S4 fits still require all declared starts, cold comparisons, full
native covariance, counts and original uncertainty distinctions. Their SEs remain
conditional on selected ecological scales. Estimated field ranges, paired GDR
spatial scores, prepared point-only search, physical state-identity tables and
integration-accuracy/recovery assessment retain their separate development gates.

### A private conditional finite-cutoff reference (original lines 908-929)

An extension reuses the existing Poisson correction and adds a conditional NB
reference. It retains the native finite-count cutoff rather than
silently substituting an infinite-count likelihood. Let Y be the observed total,
M the missed count, lambda the expected circle abundance and p the recorded
probability. The native finite-K sum is its unbounded observed likelihood times
P(M <= K-Y | Y). For Poisson abundance M is Poisson(lambda*(1-p)); with NB size r,
its conditional size is r+Y and mean is lambda*(1-p)*(r+Y)/(r+lambda*p).
This is the original truncated sum, without renormalizing the population prior.
Caller-supplied bin probabilities retain their own integration approximation.

Independent finite-count sums check this identity, including zero detections and
the exact Poisson parent. Fixed rare-availability sections show how an unbounded
conditional observed likelihood can approach a finite limit while a fixed cutoff
loses mass. Those sections shift log abundance and the removal logit intercept
in opposite directions. Cue probabilities retain their logistic law; the rare
limit does not replace that law with a log-linear hazard. Positive numerical curvature induced by that cutoff is not a finite-
intensity identification certificate. These are conditional fixed-parameter
references, not Gaussian-integrated likelihoods, derivative qualification,
zero-variance-field admission or a prepared fitter. Existing public targets,
priors, fitted parameter counts and uncertainty remain unchanged.

## From `vignettes/feature-compatibility.Rmd` (as of commit bfdc133e3)

### Bayesian, private marginal, and joint mixed admission requirements (original lines 82-89)

Bayesian preparation, targets, prediction, conditional coefficient proposals,
and covariate curves use the same native response products with frozen contrast
and fixed-transformation metadata. Private marginal/factorized components retain
their structural and numerical admission requirements. Joint mixed models use native response products and factor coding for their
supported count families, with independent cold/native/Hessian audits. Their
raw fixed-transform inputs are retained for effect curves; scale uncertainty
remains unavailable. GAMs use native refitting, preserving their smooth
and penalty targets.

### Fixed-knot splines do not qualify GDR TMB fitting or held-out scoring (original lines 139-142)

Ordinary retained training recipes continue to govern new points and raster
chunks. This response-design support does not establish GDR TMB fitting or
paired held-out scoring, which require their own likelihood and integration
qualification.

### Bayesian effect panels and gated private SQ limiting routines (original lines 346-354)

For Bayesian fits, `plot_marginal_effects(bayes_fit)` uses stored
response-coefficient draws at hypothetical predictor values under the fixed reference
normalization. Its bands are central credible intervals. It sets response
random effects to zero and does not reevaluate the landscape at every kernel
draw. For location-specific posterior expectations, use
`predict(bayes_fit, ..., kernel_inputs = new_preparation)`; use `simulate()`
when observation variation is needed. A selected-scale raster cannot replace
posterior integration over scales. Private unpenalized marginal-likelihood SQ
limiting routines remain gated pending their separate rank and boundary audits.

### One-kernel NB2 scale_inference() assessments are not calibrated tests (original lines 407-413)

One-kernel NB2 `scale_inference()` assessments retain native products and factor
contrasts in constrained, unrestricted, and bootstrap refits. The zero-effect
diagnostic removes every response column involving that kernel. It is a
feasible likelihood ceiling, not a calibrated test for an interaction or kernel
effect. Existing kernel-mean coverage studies do not establish interval coverage
for weighted SQ or interactions. All starts, failures, cold checks, dispersion
boundaries, and uncertainty qualifications remain part of the assessment.

### Checked analytic raster prediction and the normalization engine (original lines 429-434)

Checked analytic raster prediction retains every native coefficient column,
including kernel-by-factor products and fixed learned transformations. Raster
standard errors remain conditional on the selected scales. The normalization
engine uses one registered native call with stable max shifting; independent
native likelihood and cold-start checks remain in place. No comparative speedup
is claimed from validation performed alongside other R work.

### Guards and safeguards for kernel arithmetic (original lines 455-461)

Expressions must be pointwise functions of optimized kernels and literal
constants. Fixed-data inputs inside a dynamic expression, learned dynamic
bases, arbitrary or nondifferentiable functions, and kernel-valued offsets
retain explicit guards. Nonfinite domains are errors, with no clipping.
All fitting, covariance, native refit, scale-profile and private-admission
safeguards still apply. Supporting a formula does not establish identification
or calibrated intervals.

### Qualified fixed transformations and the one-case conditional-logistic adapter (original lines 489-506)

Known qualified fixed transformations use the evaluated training metadata in
profile GLM/NB2 and prepared response/scale formulas. Subset and new prediction
rows retain the original center and spread for `base::scale()`. Fixed logical
expressions inside `I()` retain R meaning during preparation; formula random
terms remain subject to the existing additive-intercept contract.

The supported one-case conditional-logistic adapter preserves original fixed
training metadata when removing the choice-set stratum. Qualified fixed
`base::scale()` and `stats::poly()` calls keep their original expression names,
contrasts, offsets, and coefficient geometry. This does not admit additional Cox
structures or provide joint scale uncertainty.

Fixed training transformations also apply to new rows with a constant predictor,
reordered sites, saved prediction snapshots, and raster prediction blocks.
Qualified `base::scale()` calls use the original training center and spread.
Native fixed spline terms retain their training knots and boundaries. These
checks preserve coefficient-conditional uncertainty; they do not establish
joint scale uncertainty or admit additional model structures.

### Metadata recovery for fits saved before the prediction-metadata repair (original lines 548-553)

Full analytic fits saved before this repair can recover prediction metadata
from their retained evaluated native frame. Numerical-only snapshots have no
such frame, so prepare them again from the full fit. This recovery does not
supply missing raw variables for a curve, change a prior, add joint scale
uncertainty or qualify a new model component. Unsupported native GAM layouts
that require a metadata change remain explicitly guarded pending component tests.

### Unmarked curve prediction contracts do not establish calibration (original lines 568-571)

Ordinary unmarked curve lists carry a `prediction_contract` recording zero
random effects, no random-effect integration and native confidence limits
conditional on the selected ecological scales. This metadata does not establish
fit eligibility or interval calibration.

### Native mgcv::bam() software checks, family coverage, the mgcv::nb() bookkeeping discrepancy, case weights, and the legacy NB family (original lines 700-849)

Package predictions from native `mgcv::bam()` fits retain the fitted factor
coding when session contrasts change. Native block controls and optimized
discrete prediction remain in place. Native parallel prediction workers also
receive saved coding and explicit term exclusions. Their mgcv version and
prediction implementation must match the parent; a mismatch stops prediction.
The checks use declared Gaussian software fixtures, not a speed benchmark or
an interval-calibration study. Mean standard errors keep the native coefficient
covariance target and omit ecological-scale uncertainty and observation noise.
Numerical snapshots cache the requested outputs without retaining a worker pool.


For ordinary BAM prediction, `type = "iterms"` returns a row-by-term matrix,
including its native intercept constant. Worker prediction retains that
geometry, selected or excluded terms, and the `iterms.type` control. These
routes also support native Gaussian BAM fits without factors and fits whose
formula contains only smooth terms. The optimized discrete provider resets
`iterms` to `terms` and issues its native warning. Use `discrete = FALSE` when
you need ordinary `iterms`; for a numerical snapshot, bind this choice with
`native_args = list(discrete = FALSE)` and `component = "iterms"`.

An `unconditional = TRUE` request does not guarantee an available
smoothing-adjusted covariance. Native BAM may discard that covariance and
warn that the correction is unavailable. Its coefficient-covariance decision
is preserved, and warning conditions from ordinary prediction workers also
reach the parent call. Read these SEs as conditional on the ecological scales;
they do not include observation noise or full-posterior uncertainty. Native
EDF and likelihood parameter counts remain distinct from the two selected
ecological parameters in the declared Gaussian software checks.


Native mgcv versions may also warn when selecting or excluding terms from a
model with no parametric coefficients. That provider warning is retained
alongside the native numerical result. It is separate from the fitted-mode
convergence, score and curvature checks.


Warnings transported from ordinary BAM workers retain their native messages,
classes and prediction controls. Their calls identify the native method and
refer to fitted models and prediction rows symbolically, avoiding extra copies
of those data in warning payloads. Provider version and the symbolic argument
names remain attached as diagnostic context. The worker's original warning
is retained; unrelated native warning calls are unchanged.


Native Gaussian, binomial, Poisson, Gamma and negative-binomial BAM fits use the shared point,
numerical-snapshot, checked-map and raw-reference curve routes. Ordinary
prediction supports factor interactions, `I(z^2)` terms and native smooths;
set `discrete = FALSE` when requesting ordinary `iterms` output. Native
parallel scalar and matrix outputs retain their provider geometry, which can
differ from serial geometry. Term matrices retain their labels and intercept
constant. Native effective degrees of freedom supply the coefficient and
smoothing contribution to AIC/BIC, with the actual family nuisance parameters
described below; selected ecological scales and applicable expow shapes are
added once.

For a single-predictor GAM or BAM, `plot_marginal_effects()` inherits the GLM
automatic response transformation. Both values of `link` retain response
curves, with inverse-link bounds from conditional link-scale Wald intervals.
Use `predict(..., type = "link")` for link-scale predictions. Multi-component
GAM curves have an explicit component choice: `link = FALSE` selects a linear
predictor and `link = TRUE` selects a native response component. These
conventions preserve the existing interfaces.

Binomial response predictions are probabilities, including models fitted to
success/failure counts. Trial counts enter the fitted likelihood and
covariance; they do not automatically convert predictions to counts for new
trials. Poisson response predictions are conditional expected counts. Neither
response SEs nor these curve intervals add ecological-scale uncertainty,
observation noise or full-posterior uncertainty. Native BAM rejects general
families; use an appropriate native GAM for multi-predictor families instead.


Optimized discrete prediction (`discrete = TRUE`) is also checked for native
binomial and Poisson BAMs across all four kernels, with point, saved-snapshot
and raster output compared to the same native prediction batch. Discrete
predictions can depend on the batch used to construct their approximation;
ordinary predictions or a different row batch are separate targets. Native
warnings, term matrices and conditional coefficient uncertainty are retained.
These software checks do not establish statistical calibration or speedups.


#### Native BAM family coverage

The native-refit adapters also support Gamma models with a log link and
negative-binomial models using `mgcv::nb()` with fixed or estimated `theta`. The software checks
pair a kernel mean with weighted SQ, fixed factor coding, a numeric product,
a polynomial term and a native smooth. They retain each fit's convergence,
coefficient-score, information and covariance checks.

| Family | Ordinary prediction | Optimized discrete prediction | Information-criterion nuisance parameters |
|:-------|:--------------------|:------------------------------|:-----------------------------------------|
| Gaussian | Points, snapshots, maps, curves and workers | Native controls and checked outputs | Native effective degrees of freedom and estimated residual scale |
| Binomial or Poisson | Points, snapshots, maps, curves and workers | Same-batch points, snapshots and maps | Native effective degrees of freedom; no residual scale |
| Gamma with log link | Points, snapshots, maps, curves and workers | Same-batch points, snapshots and maps | Native effective degrees of freedom and estimated residual scale |
| Negative binomial | Points, snapshots, maps, curves and workers | Same-batch points, snapshots and maps | Native effective degrees of freedom and estimated `theta`; fixed `theta` adds none |

The Gamma and estimated-theta negative-binomial cases cover Gaussian,
exponential, exponential-power and fixed kernels. Fixed-theta cases use a
Gaussian kernel. These tests check adapter compatibility at declared spatial
scales; they do not demonstrate ecological-scale recovery or interval coverage.
Some ordinary family fixtures require a tighter native fitting tolerance to
meet the unchanged coefficient-score gate. Package fitting defaults remain
unchanged; inspect diagnostics for your own fits.

Gamma response predictions are conditional positive means. Negative-binomial
response predictions are conditional expected counts. Their mean SEs and curve
bands use native coefficient covariance, holding ecological scales fixed.
They do not add observation noise, uncertainty from estimating `theta`,
ecological-scale uncertainty or full-posterior uncertainty.

Some native mgcv versions label `mgcv::nb()` BAM scale as estimated even
though the fitted family has unit scale. In those verified cases, native
`logLik()` and its parameter count each include one extra unit. `aic_tab()`
and `bic_tab()` check the normalized negative-binomial response density and
the actual native effective degrees of freedom before correcting that specific
bookkeeping discrepancy. Estimated `theta` still counts once; fixed `theta`
does not count. Other discrepancies produce an error. The correction leaves
AIC unchanged, while AICc and BIC use the corrected count. Fitted coefficients,
covariance and native diagnostic objects are retained.

An affected table has a `likelihood_adjustments` attribute recording the
original input position, native and used log-likelihoods, regression parameter
counts, and nuisance counts. These recorded regression counts precede any
selected ecological scales and applicable shapes added to the table's `K`.
A fixed radius contributes no fitted parameter. A selected hard cutoff still
counts once even though its Hessian radius SE is unavailable.


#### Case weights and the legacy NB family

Native BAM case-weight checks include zero and fractional prior weights.
The fitted covariance and native observation count are retained. Zero-weight
rows remain part of the declared preparation and its frozen reference
normalization; they do not supply likelihood information. Prediction SEs still
condition on selected ecological scales and native nuisance estimates.

For the legacy `mgcv::negbin(theta)` family, supply `scale = 1` explicitly when
you want the checked known-scale NB target. The BAM default may estimate a
working scale for this family. That default is a different native convention;
the adapter does not silently replace it. Known fixed theta and known unit
scale add no fitted parameters. The legacy provider keeps its native dispatch.

With very large weights, subtracting two large log-likelihoods can obscure a
one-unit bookkeeping discrepancy. Verification then checks the original native
AIC arithmetic and compares rescaled response densities at the same tolerance.
Rescaling is used only for verification. Reported likelihoods use the original
weights, with the same parameter counts and covariance. An affected adjustment
records `verification = "native_aic_and_rescaled_density"`. Altered metadata
continues to produce an error. Arithmetic copies used to test this rounding
case are not refitted models or evidence of inferential calibration.

### Explicit-fold cross-validation qualification statements (original lines 851-896)

#### Explicit-fold GLM validation

`cross_validate()` provides training-only scale refitting and held-out scoring
for Gaussian identity, Poisson/NB2 log, binomial logit and Gamma log GLMs.
Gamma and fixed-theta NB2 require native `refit`; estimated NB2 supports
`profile_glm`, `joint_nb2` and `refit`. It accepts explicit
analysis/assessment row IDs, retains failures and uses frozen training bins and
normalization. Cell-level inputs are required; lean-only preparation is rejected.
Modern grouped raster inputs retain validated per-row fixed-map identities in
each split; they do not introduce annual random effects or varying scales. Supported kernels and metrics follow the selected
backend's existing guards. Training-only theta/Gamma dispersion and normalized
log scores preserve fixed/estimated nuisance counts and native coefficient SE
conventions. NB2/Poisson count comparisons and Gamma continuous comparisons
retain likelihood-measure guards. Selected folds record their actual family and
nuisance parameters. This GLM route does not add
GAM, mixed-model, unmarked or Bayesian scoring. The separate conditional-choice
route is described below. See
`vignette("inference-contracts", package = "multiScaleR")` for the workflow.

`grouped_cv_folds()` and `spatial_cv_folds()` generate compatible whole-unit
splits, with optional spatial training buffers and shared-group block joining.
`cross_validate_glm()` supplies ordinary GLM reference models. `cv_tab()` ranks
matching held-out scores; `select_cv()` evaluates training-only AIC/AICc/BIC
selection with retained fold fits and original parameter counts. These helpers
do not expand the supported scoring model classes or any private admission.
Worked examples are in `vignette("cross-validation", package = "multiScaleR")`.

For memory-efficient training selection, request `keep_selection = TRUE` in
each CV call. Compact records retain native likelihoods, complete parameter
counts, training identities, controls and diagnostics. Full models remain
optional via `keep_fits = TRUE`. Record-based selection retains the existing
comparability, failure and uncertainty contracts; it does not qualify additional
model classes for CV scoring.

`cross_validate_clogit()` qualifies a separate conditional-choice CV target for
complete one-used-step sets, with native references or supported scalar spatial
backends. Every set stays wholly in training, assessment or exclusions. Scores
average natural-log choice loss over sets while retaining all alternatives.
Native exact/Efron/Breslow methods, event-based nobs and complete K remain in
training selection, including compact choice records. Generic censored Cox,
case weights, robust/penalized terms and other model-class scoring still need
their own qualification. The CV-only null likelihood adapter checks survival's
stored scalar against the complete-set offset likelihood; it does not change
other native Cox objects. Priors, numerical gates, private admission, conditional
SEs and uncertainty distinctions remain unchanged. See
`vignette("cross-validation", package = "multiScaleR")` for a worked example.

## From `vignettes/cross-validation.Rmd` (as of commit bfdc133e3)

Development-status passages removed from the vignette during the rewrite, verbatim.

### Scoring adapters that still need qualification (original lines 64-67)

The GLM scoring adapters support Gaussian, Poisson, binomial, NB2 and Gamma
targets. Complete matched choice sets use the separate conditional-choice
adapter described below. Mixed, GAM, unmarked and Bayesian scoring require
their own qualified predictive targets and adapters.

### Prediction repair for learned transformations and its limits (original lines 126-140)

For spatial CV, these constants use the final training predictors at the fitted
ecological scale. Initial-scale predictors and assessment buffers cannot set
them. Native formulas, coefficient estimates and SEs, likelihoods, parameter
counts, optimizer diagnostics and uncertainty targets stay unchanged. This
prediction repair also applies to retained fold fits and complete-set choice CV.
It does not widen a backend's admitted spatial expressions or fitting methods.

Automatic recovery is limited to these known learners with numeric inputs and
supported pointwise arguments, such as `scale(log(x))`. It does not infer how a
custom function learns from data or re-execute arbitrary functions to obtain
metadata. For other learned operations, including a formula containing
`mean(x)`, train an explicit preprocessing recipe within each analysis split
and apply its saved constants to the assessment rows. Unsupported arguments to
a known learner leave a failed fold with its original error. Check `complete`
and `errors` before comparing candidates.

### Joint-NB2 boundary gate and initializer (original lines 222-229)

Native convergence or theta-estimation warnings retain failed folds; do not
drop them to obtain a ranking. The exact Poisson limit follows the original
`joint_nb2` boundary gate and flag, rather than replacing theta with a large
finite number.
Joint-NB2 initializes coefficients with a native fixed-size GLM at `theta = 1`,
then estimates dispersion with the existing joint engine and its boundary
checks. This avoids requiring finite-theta alternation before the exact boundary
can be assessed. The initializer's size is neither a fixed fitted theta nor a prior.

### Known Gamma dispersion and native AIC (original lines 253-258)

Use `theta = 3` or `dispersion = 0.4` only if that parameter is known or fixed
independently of assessment outcomes. Known Gamma dispersion uses its declared
normalized density from the start of fitting. Merely setting a native Gamma
family's dispersion field would change its parameter count while leaving its
AIC function's deviance estimate in place. The CV adapter avoids that mismatch
without changing global GLM methods or existing fitted objects.

### NB2 and Gamma comparison guards and backend admission (original lines 260-275)

Both families require integer nonnegative frequency weights. Native nobs is
the number of positive-weight training rows, not the sum of their frequencies.
Gamma outcomes must remain strictly positive and NB2 counts nonnegative
integers, including zero-weight rows. Their default scores are RMSE, MAE and
negative log score. Deviance is not offered as a comparison loss because its
saturated-density constants change with the nuisance parameter. RMSE/MAE can
compare compatible conditional means across families; density scores and
information criteria require the same continuous or discrete observation
measure. Gaussian density scoring is not qualified here.

For spatial fits, estimated NB2 supports `profile_glm`, `joint_nb2` and `refit`.
Gamma and fixed-theta NB2 require `refit`, which is selected when their backend
argument is omitted. Explicit incompatible choices stop before fitting.
Existing kernel/metric and numerical guards remain in force. Mean and supported
weighted-SQ predictors use the same frozen training normalization and reference
contract as the other GLMs.

### Scope of plug-in scores (original lines 287-289)

These scores condition on fitted mean and nuisance parameters. They integrate
neither parameter nor selection uncertainty, provide no observation intervals
and establish no calibration claim. Existing priors and diagnostics are unchanged.

### Test-suite coverage of training-only scale refitting (original lines 323-327)

Reuse the same explicit folds with raw-cell preparation. This workflow is not
run during vignette builds because it adds a scale optimization for every fold.
The test suite checks training-only scale refitting and score agreement with
native GLMs. Choose the extraction domain and initial scales independently of
assessment outcomes.

### Raw-unit predictors, backend restriction and unchanged gates (original lines 383-388)

Raw mode currently requires the native `refit` backend. Analytic backends retain
their complete candidate-standardization derivatives, and Bayesian methods keep
their frozen reference targets. Parameter penalties, retained fold failures,
conditional coefficient uncertainty and the absence of selection/parameter
integration in these scores are unchanged. This option does not add mixed-model
GDR held-out integration or new uncertainty intervals.

### Cell geometry and the legacy strict gate (original section "Retain cells intersecting the extraction buffer", lines 418-428)

Native raster preparation records both cell dimensions. A cell can intersect
the extraction buffer even when its center lies beyond the nominal radius.
CV retains those values and fixes the bin domain using the requested radius
plus half the pixel diagonal, rather than a maximum learned from held-out rows.
Representative distances and normalization are still learned only from training
rows. Geometry is frozen with the bin reference. Legacy saved inputs without
cell geometry keep the original strict gate; prepare them again with the current
constructor if that gate rejects their cell-center distances. Extraction,
kernel/metric definitions and fitted parameter bounds are unchanged.

### Conditional-choice null likelihood validation (original lines 475-481)

The null model has no estimated response coefficients. CV validates survival's
stored scalar null likelihood against the complete-set offset likelihood before
using it for training selection. Missing or inconsistent metadata remains a
failure. Other native Cox objects and methods are unchanged. The native training
table uses one event per analysis set for AICc/BIC and counts all estimated
response coefficients and ecological scale/shape parameters. It does not use
the number of available alternatives as the event count.

### Conditional-choice spatial backends and adapter limits (original lines 491-498)

To estimate ecological scales within those same folds, supply raw-cell
`kernel_inputs` and choose `backend = "conditional_logit"` or `"refit"`.
Existing kernel/metric admission, independent starts, cold/native checks and
curvature rules remain. Initial and fitted normalization and representative
bin distances use analysis rows only. Assessment predictors use those frozen
references, including when their distances enter different training bins.
Prepared raster-map groups remain fixed covariate inputs and add no annual or
random effects by themselves.

### Conditional-choice adapter limits (original lines 507-514)

Use the same declared exact, Efron or Breslow method across candidates. This
adapter accepts unweighted, unpenalized one-case sets. Robust clustering,
frailty/special Cox terms and generic censored survival prediction need their
own scoring targets. `cv_tab()` keeps conditional choice separate from GLM
scores. As with the GLM workflow, inspect original diagnostics, retain failures
and use outer assessment if these results tune the candidate set or selection
rule. These point scores integrate no parameter/selection uncertainty and
provide no observation intervals or calibration certificate.

### Gaussian response-field CV, uncertainty scope and start gates (original lines 554-575)

Field states are integrated jointly at the fitted training parameters. Parameter
uncertainty, including uncertainty in correlation range, observation SD and any
selected ecological scale, is not integrated into this plug-in density. Native
mean-prediction SEs keep their own meaning and are not observation predictive
variances. New-location field variance, remaining training-field uncertainty
and observation noise are identified separately in the returned predictions.
Software agreement is not a calibration or recovery assessment.

Supplying `kernel_inputs` also re-estimates its ecological scales inside each
training fold. Raw cells, bins, normalizers and supported scale/poly/ns/bs
recipes use analysis rows only. Assessment and excluded buffer responses never
enter fitting. The response-field range and ecological kernel bandwidth remain
separate estimated quantities. Both native covariance starts must pass the
original convergence, score, curvature and agreement gates. Failed folds and
their compact audit records are retained; no passing replacement start is used.

Compact training records always support AIC/AICc/BIC selection for these field
results. They count native regression, range and variance parameters plus
estimated ecological scales, while integrated field realizations add no penalty.
Choosing a model from held-out ranks still needs independent outer assessment.
The adapter estimates no Matern shape, additional random effects or paired GDR
field score. Those routes retain separate qualification guards.

### Selection records and metadata checks (original lines 588-592)

`select_cv()` checks those records against the split and target identities and
uses the same likelihood-comparability guards and criterion table constructors.
If you retain both records and full fits, they must agree exactly. A metadata
failure is retained as a selection failure and does not erase a candidate's
otherwise valid held-out predictions. Records do not change diagnostic eligibility.

## From `vignettes/distance-removal-audits.Rmd` (as of commit bfdc133e3)

The vignette was rewritten as a task guide for fitting distance-removal models
with scales of effect. It keeps one "Current restrictions" table. The passages
below are quoted verbatim, including their code chunks. Three functions named
here (`audit_gdr_integration()`, `gdr_start_values()`, `score_simulation()`)
are scheduled to be un-exported and no longer appear in the vignette.

### Scoring simulated estimates with `score_simulation()` (original lines 127-155)

A density coefficient estimated on the log scale must be exponentiated before
it is compared with response-scale density truth. Transform interval endpoints
as well. Coefficient positions are not a safe way to select targets when model
formulas have different numbers of terms.

```{r score}
attempts <- data.frame(
  replicate = 1:4, target = "density", units = "birds/km2", truth = 200,
  estimate = c(log(180), log(200), log(220), NA),
  lower = c(log(150), log(170), log(185), NA),
  upper = c(log(215), log(235), log(260), NA),
  eligible = c(TRUE, TRUE, TRUE, FALSE),
  reason = c("", "", "", "native_hessian_failure")
)
scores <- score_simulation(attempts, input_scale = "log", level = 0.95)
scores$summary
scores$attempts[c("replicate", "score_eligible", "score_reason", "interval_status")]
```

All four attempts remain visible. Bias, RMSE, and variation condition on the three
eligible point estimates. Interval coverage uses eligible finite intervals and
reports its own denominator. These intervals are supplied examples; the helper
does not construct or validate a 95% procedure. Retain original convergence,
Hessian, boundary, and other diagnostic gates when defining eligibility.

Relative bias is undefined at zero truth. CV is undefined at zero mean or with
fewer than two eligible attempts. Set `identified = FALSE` when scoring ecological
scales under a zero habitat effect: an optimizer can return a scale even though
the generating effect does not identify one.

### Response definitions and individuals that never produce a record (original lines 206-217)

Paired first cells, native two margins and complete linked histories are
different responses. Their likelihood constants and conditioning differ.
The helper does not authorize comparing their objectives in one AIC table.
Five within-count bands also do not constitute five independent primary visits.

Better description of recorded histories does not by itself identify birds
that never produce an eligible record. Multiplying all observable probabilities
by `1 - w` and expected abundance by `1 / (1 - w)` preserves the observed
Poisson or NB thinning distribution. A positive Hessian or finite optimizer
bound cannot resolve that exact equivalence. Keep target-population assumptions,
observation-process checks and any external calibration explicit when reporting
density; conditional recorded-count prediction is a different validation target.

### What interval preparation does not establish; the first-interval law (original lines 274-289)

The index is retained as an observation predictor, so a justified native
removal formula can use `~ I(interval == 1)`, `~ factor(interval)`, or an
interaction with an observer covariate. These examples illustrate data layout;
two events do not establish that a removal model is identifiable. Native fitting
and `aic_tab()` retain the parameters implied by the actual component formulas.
Preparation does not qualify an ignored offset, conditional-random prediction,
another likelihood, or a density correction.

For equal unit intervals, an initial encounter probability pi followed by
activation probability b gives the simple Markov first-time law regardless of
persistence after an encounter. The existing native first-band interface can
represent that law. First times alone cannot estimate persistence or justify a
clock-response interpretation. Use linked complete histories when studying
re-detections, and qualify assumptions about unrecorded individuals before
interpreting corrected density. Unequal native durations retain per-unit
probability meanings; this helper labels their units without converting them.

### Auditing integration at the refitted coefficients (original lines 294-344)

An original fit can pass a latent-count cutoff check while a resampled or
scale-refitted fit fails it. `audit_gdr_integration()` evaluates the stored
coefficients without optimizing them. For a qualified single-primary native
C Poisson or NB fit, request an explicit cutoff ladder and inspect the original
cutoff alongside the additional evaluations.

```{r native-integration-audit, eval=FALSE}
# native_fit is the fitted model being assessed, with its original data/design.
integration <- audit_gdr_integration(native_fit, cutoffs = c(80, 160, 400))
integration$stored_point
integration$cutoffs
integration$references
integration$distance
# Additional exact finite-cutoff reference for Poisson abundance.
integration$finite_poisson
```

`truncation_loss` compares each finite-cutoff negative log likelihood with an
infinite-sum observed law using the same native distance bins. A failing loss
is a numerical integration finding, even when local convergence and Hessian
checks look satisfactory. The helper records the requested evaluations and
their failures. It leaves the fitted cutoff and all coefficients unchanged;
it does not select a replacement fit or change an optimization target.

For Poisson abundance, `finite_poisson` also evaluates an exact finite-cutoff
identity with the same native distance/removal probabilities. The observed total
has mean `lambda*q`, with a Poisson CDF correction for unobserved individuals up
to `K-n`. This reference removes its latent-abundance sum while retaining the
original finite cutoff and multinomial constants. It is not a renormalized
truncated abundance distribution. The native cutoff ladder still runs, and a
passing identity comparison does not replace a failing truncation-loss check.
Inspect `comparison$agrees` and each retained error/warning record. NB is
explicitly unqualified for this additional identity; its original audit remains.
A longer cutoff ladder retains compact summaries instead of a set of event-level
predictions for every cutoff. This reference supplies neither a TMB refitter
nor Gaussian integration or SEs.

The separate distance comparison uses stable analytic annulus integrals.
Native unmarked 1.5.2 exponential bins use 100 trapezoids per bin, while native
half-normal arithmetic can lose precision near the uniform limit. Differences
between those bins and the stable analytic reference are distinct from finite
latent-count truncation. Invalid native probabilities remain failed records.
Increasing the cutoff cannot repair a different distance-integration target.

This audit covers complete single-primary fixed-effect C models with Poisson
or NB abundance and half-normal, exponential or uniform point-transect keys.
Other engines, random effects, missing histories and multiple primaries need
separate qualification. Numerical agreement does not establish population
density identification, biological fit or valid uncertainty. Native coefficient
intervals still condition on the selected ecological scales, and all parameter
penalties and existing fitting guards remain in place.

### Native fitting and remaining integration gates (original lines 349-403)

Use `unmarked::gdistremoval()` for native fitting and its component-specific
prediction methods. Data and recovery checks do not replace that likelihood or
admit a new fitting backend. Check the provider's likelihood before using a fit
for ecological-scale optimization, information criteria, or mapping.

multiScaleR rejects two unqualified GDR combinations. In an installed unmarked
1.5.2 fixed-point audit, TMB ZIP likelihoods ignored the zero-inflation coefficient
and instead used probability 0.5. The C engine's abundance, distance, and removal
offsets did not change its likelihood. Active availability offsets are also held
pending provider qualification. These restrictions apply to scale refitting,
AIC/BIC, spatial-data preparation, prediction snapshots, maps, and effect curves.
An availability formula is inactive for one primary period; a distance formula is
inactive for a uniform distance key. Actual R offset terms are distinguished
from ordinary expressions such as `I(offset(x))`.

An offset-free C ZIP model keeps its native likelihood and all fitted parameters.
It is not a fallback for a random-effects model, which requires the TMB engine.
Do not drop an intended offset, change its penalty, or switch engines simply to
bypass a restriction. A different provider version needs new qualification.
The native ZIP ancillary back-transformation in the audited version returned
odds rather than probability, so its native ZIP simulation path is not qualified
here. `sim_dat_gdr()` remains restricted to Poisson observations.

Native TMB fits request `sdreport()` even when a caller requests `se = FALSE`.
That setting therefore does not establish a faster Hessian-free TMB refit.
Preserve final covariance diagnostics and distinguish native optimizer starts
from the ecological-scale optimizer's starts.
 Native optimizer convergence,
finite Hessian uncertainty, and robustness to declared starts are separate checks.
Keep the same survey definitions and native integration cutoff K when comparing
starts within an ecological-scale optimization. Study K sensitivity separately.

GDR information criteria compare both response histories and the stored survey
protocol, including primary-period layout, distance bins and units, and removal
durations. Equal distance/removal totals do not establish comparable data.
Use the explicit event audit before comparing fits. Native frames do not retain
a removal-time unit label, so matching stored numbers alone cannot verify that
one analysis used minutes and another used seconds. Keep units and event
definitions consistent, and assess integration cutoff K separately.

Only active component formulas identify ecological spatial scales. A one-primary
model does not estimate availability, and a uniform key does not estimate a
distance scale. Variables occurring only in those inactive formulas are excluded
from spatial-variable discovery. Optimized static site predictors also need
distinct names from yearly and observation predictors. This prevents a repeated
native table from shadowing the site values being updated. Dot formulas expand
against the original active component's training-column names; transformations,
factor names and interaction inputs are identified through R's formula terms.

For density maps, specify density units and the surveyed area, and distinguish
expected ecological density from observed counts. In a zero-inflated model,
the count-component intensity differs from the unconditional expected density.
Conditional coefficient uncertainty also differs from uncertainty that includes
estimated ecological scales. Neither a profile curve nor a coefficient-only
interval automatically provides scale-inclusive coverage.

### Retaining numeric constructor starts (original lines 407-428)

`gdr_start_values()` reads the stored fixed optimizer values and original
constructor order. Native TMB input order can differ from fitted order, including
the position of dispersion and log-SD blocks. The returned map identifies the
component and within-component position of every coordinate; repeated labels
are not used to guess a match. Latent random modes stay separate.

```{r native-start-metadata, eval=FALSE}
# native_fit is an existing user-owned unmarkedFitGDR.
start_record <- gdr_start_values(native_fit)
start_record$map
constructor_starts <- start_record$starts
```

The plain numeric/list record can be saved without a model or executable tape.
Its fingerprint describes coordinate metadata, not the full data, likelihood
or prediction context. Restore an objective only from its original data,
transformations, formulas and evaluated controls, with all original checks.
The stored convergence code and objective retain their original meaning; finite
starts do not establish convergence, covariance or identification. The coordinate
count does not replace the native AIC penalty. No likelihood, score, Hessian or
covariance is evaluated by this inspector. In the inspected native TMB path,
top-level controls and `se = FALSE` do not provide a likelihood-only fit.

### The native refitter contract and TMB settings (original lines 430-490)

An opt-in native C or TMB refitter is available through `make_gdr_refitter()`.
The C path recovers its original evaluated settings. The TMB path requires
explicit evaluated settings and keeps the fitted native engine. It supports
native Poisson/NB targets and normalized Gaussian random effects subject to
the existing offset, zero-inflation and covariance guards.
Its `contract.rds` binds numeric starts to the original data/design fingerprints,
evaluated settings, standard-error request and native provider identity. Every
refit records the actual permitted input design and a contract-linked input
fingerprint. These hashes identify the retained inputs; supplied provenance
still needs independent verification.

Each native attempt writes a request before entering the constructor and a
separate completion afterward. A request without a completion documents an
unfinished attempt. Retain it. Completed records include the optimizer's numeric
coordinates, objective, evaluation counts, convergence code and message when
available, including failed diagnostics. The records do not contain raw data,
design matrices, Hessians or executable objectives. A saved record cannot replace
fresh native construction, independent starts or the original final uncertainty
checks. No new coefficient score is evaluated by checkpointing.

Create it for the original native model and keep its audit directory. C cold
starts retain their named-coefficient interface, skip intermediate coefficient
Hessians and restore the original SE request. TMB starts include every fixed
coefficient and log-SD in native constructor order, excluding latent modes.
Do not substitute the fitted optimizer's block order or repeated coefficient
labels for the constructor map. Each TMB attempt starts its latent coordinates
at zero and creates a fresh native objective.

```{r tmb-refitter-settings, eval=FALSE}
# Inspect the native constructor coordinates before changing a cold start.
coordinates <- gdr_start_values(native_tmb_fit)
cold <- coordinates$starts
cold[coordinates$map$kind == "log_sd"] <- log(0.4)
# Also change an appropriate fixed coefficient if the model has no log-SD.
cold[1L] <- cold[1L] - 0.2
refit <- make_gdr_refitter(native_tmb_fit, list(cold = cold),
  audit_dir = tempfile("gdr-tmb-audit-"),
  tmb_settings = list(method = "BFGS", se = TRUE, threads = 1L,
    control = list(maxit = 2000L, reltol = 1e-10)))
# Supply refit as refit_fn to multiScale_optim(backend = "refit").
```

Declare the settings deliberately. A saved call does not establish the original
evaluated TMB controls, and unmarked 1.5.2 drops optimizer extras on its TMB path.
The adapter forwards the declared controls locally without patching unmarked.
It evaluates native `sdreport()` at the returned optimizer coordinates so that
the native estimates, covariance and optimizer metadata refer to the same point.
It retains finite marginal scores, native positive-Hessian status and full
joint-precision checks. Bad covariance remains a failed fit.

TMB still computes native `sdreport()` for intermediate fits, even with
`se = FALSE`. This interface supplies working native refits, not a claim of
Hessian-free TMB optimization. Final covariance keeps fixed/latent and component
cross-covariances; ecological-scale and model-selection uncertainty are not
added to it. Each request, completed attempt and final outcome
has an atomic compact RDS record. Higher likelihoods with failed convergence or
unusable final covariance remain failures rather than being replaced silently.
Input fingerprints are computed in memory; raw data and candidate-model arrays
are not checkpointed. These native coefficient intervals remain conditional on
the selected ecological scales. Native mixture codes also distinguish default
AIC/BIC labels without changing their parameter penalties.

### Raw covariate units, bound-limited scales and remaining gates (original lines 492-524)

For an initial native model whose optimized habitat predictors are in their
kernel or metric units, `multiScale_optim(..., backend = "refit",
covariate_transform = "raw")` preserves those units throughout refitting and
spatial prediction. The default continues to standardize each candidate. For
example, a knot at habitat fraction `0.2` has a different meaning from a knot at
`0.2` candidate standard deviations. Fit the initial model and specify its
formula in the intended units; the optimizer does not reinterpret an existing
formula or its knots.

Raw mode stores an identity transform for optimized predictors, so a later
prediction batch cannot supply new means or SDs. It leaves upstream weighted-SQ
reference normalization, fixed-metric transforms, native observation roles,
sample identity and parameter penalties under their existing contracts.
Coefficient SEs remain native and conditional on the selected ecological
scales. Changing units does not establish learned GDR prediction contexts or
paired-process CV. Native TMB refitting has the separate contract above. Requested and actual spline knots
can differ under the native basis rules; inspect that distinction rather than
treating a near-zero knot as a biological threshold.

A smooth ecological-scale estimate at a search bound does not support an
ordinary outer Wald interval, even if its Hessian is positive. The refit backend
stores the bound and original curvature diagnostics and withholds its joint
smooth-kernel Wald block. Coefficient SEs still condition on the selected scales,
and parameter penalties are retained. Explicit profiles can describe likelihood
support, but a bound-limited profile does not resolve the missing tail. Native
prediction contracts retain both bound and extraction-window flags.

Complete GDR component prediction contexts,
paired-process cross-validation, and scale-inclusive bootstrap intervals remain
separate integration gates. Existing GLM CV and NB2 scale-inference restrictions
remain in force. See `vignette("cross-validation", package = "multiScaleR")` for
the currently supported CV targets, and `vignette("inference-contracts",
package = "multiScaleR")` for uncertainty distinctions.

### Inspecting the stored state basis (original lines 540-576)

For the existing single-primary fixed-effect native C Poisson/NB audit domain,
set `design = TRUE` to inspect the actual stored component matrices. The added
summaries report column rank, condition number and the original matrix-column
scale. They do not change fitting tolerances or numerical admission gates.

```{r gdr-design-audit, eval=FALSE}
# This illustration assumes hab is the only state predictor. Include all
# required state covariates, including factor/year levels, for your own model.
# Supply scenarios in the predictor units retained by the fitted native model.
audit <- audit_gdr_integration(model, design = TRUE,
  design_data = data.frame(hab = c(0, 0.2, 0.8)))
audit$design$components$lambda
audit$design$state$splines
audit$design$projection
```

State spline metadata are reported only after reconstruction exactly reproduces
the original native matrix, including factor interactions. Requested spline
calls and retained numeric knots are separate fields. A native constructor can
move a requested knot that coincides with a boundary, especially with zero-heavy
cover. The audit retains that distinction, its warnings and any reconstruction
failure. It never repairs knots or interprets a near-zero knot as a threshold.
Natural-spline basis columns have linear tails; added products or transformations
can change the full predictor's tail behavior. Rank and support do not establish
monotonicity or a single response peak.

Optional scenarios use the verified training prediction calls and original
factor levels/contrasts. Counts and knots use modeled predictor units; a legacy
standardized column does not automatically recover raw cover. For bare numeric
spline inputs, the audit reports counts
outside observed training ranges and outside spline boundaries separately.
Other expressions retain native matrix/recipe diagnostics with an explicit
support status. The returned summaries contain no full matrices or scenario
rows. Scenario-matrix rank describes that batch and does not reassess fitted-model
identifiability. They provide no response predictions, coefficient/kernel uncertainty,
learned GDR component context or new fitting admission. TMB, random-effect,
multiple-primary and other existing audit guards remain.

### Retaining component training recipes (original lines 580-630)

Predictions use the component's actual training rows. This matters when a
formula contains learned spline knots, an orthogonal polynomial, nested
`scale()` calls or factor interactions. Rebuilding those terms from a single
site or a raster chunk can change the fitted design. The checked GDR adapter
reconstructs the original design first, then retains its prediction recipe.
Native C and TMB component covariance and confidence limits are preserved.

| Type | Native component | Original training rows |
|:--|:--|:--|
| `lambda` | Abundance or density parameter under the fitted model | Sites |
| `phi` | Availability, when there is more than one primary period | Site by primary period |
| `dist` | Distance-detection parameter for the fitted key | Site by primary period |
| `rem` | Native removal parameter | Site by primary period by removal interval |

These outputs retain the provider's parameter interpretation. They are not
automatically survey averages, observed detection probabilities or expectations
integrated over random effects. In particular, native C ZIP component semantics
are retained. Inactive availability or uniform-distance components do not gain
an estimated response curve.

```{r gdr-trained-predictions, eval=FALSE}
# Use the covariates in the fitted model's units. Do not recenter this batch.
site_rows <- unmarked::siteCovs(msr_fit$opt_mod@data)[c(2, 9), , drop = FALSE]
# A numerical snapshot keeps these already prepared rows and controls.
snapshot <- prepare_prediction(msr_fit, newdata = site_rows,
  component = "lambda", native_args = list(re.form = NA))
predict(snapshot)

# Select an active component explicitly. Its own original rows define the grid.
plot_marginal_effects(msr_fit, type = "dist", ylab = "Distance parameter")
```

Explicit data-frame predictions prepare only the requested component. Spatial
`predict()` calls still require newdata and matching kernel_inputs
together. Use the snapshot route above for already computed native columns.
Native GDR frames retain the original site/period/interval layout and missing-input
checks. Raster chunks use the same terms and the component's original logical
and factor types. The checked original covariates are reused across chunks,
avoiding repeated reconstruction of every training observation table.

Keep the original covariates, fitted design and recoverable training metadata.
A mismatched design, unavailable factor coding or unsupported native layout
stops. Conditional learned random-design transformations remain guarded.
`re.form = NA` sets effects to zero, while ordinary conditional predictions use
the fitted known-group effects. Neither operation integrates new groups. Native
coefficient uncertainty remains conditional on selected ecological scales;
kernel/shape selection uncertainty and observation noise are not added.

This prediction adapter is separate from the numeric likelihood/design audit.
It does not enable paired-process CV, a spatial response field or the private
RTMB scale likelihood. Those targets retain their own qualification gates.

## From `vignettes/bayesian-scales.Rmd` (as of commit bfdc133e3)

### Simulation evidence behind scale formulas (original lines 186-191)

NUTS and slice sampling can fit multiple scale-group SDs. The marginal approximation currently
integrates exactly one scale-group SD and rejects models with more than one.
Additional response random intercepts are distinct from this limit. Broad
simulation evidence is strongest for one landscape predictor with annual
pooling and Gaussian kernels; a valid formula is not itself evidence of good
identification or interval calibration for a more complex model.

### Weighted SQ routes and private limiting routines (original lines 245-251)

The prepare, sampling, joint Laplace, and marginal Laplace routes use the same
SQ posterior and existing diagnostics. Numerical derivative and fit checks do
not establish posterior approximation accuracy or interval calibration. Other
optimized surface metrics and exponential-power SQ remain unsupported. Fixed-radius
metrics can enter the response or scale formulas as verified numeric columns,
retaining their original preparation means and sample SDs. The separate
private marginal-likelihood limiting routines still reject SQ.

### Features that need separate likelihood and validation work (original lines 285-292)

The landscape coefficient currently stays common across years even when its
scale varies. Kernel-by-year interactions and response random slopes are not
supported. The interface also lacks correlated random effects, temporal
autoregression, spatial random fields, zero inflation, custom families and
joint detection or multispecies likelihoods. These would require explicit
likelihood, design and validation work. Bayesian inference alone does not add
them to the package. Existing likelihood wrappers support some models outside
this Bayesian interface, so it is not a replacement for every existing backend.

### Hierarchical scale models as a development capability (original lines 301-302)

Hierarchical scale models are new capabilities of this development interface;
they are not exclusive to Bayesian inference in statistical theory.

### Retries of the stored example fits (original lines 516-519)

These checks use all retained chain iterations, rather than treating the
smaller bundled draw count as additional evidence. If a longer retry was
needed, its actual iteration count appears in the table and the original
failed attempt remains in the development study.

### Statistics used by the comparison studies (original lines 756-758)

The comparison studies use log-distance errors and distance quantiles, so their
primary recovery and calibration statistics do not require a finite distance
variance. These clarifications leave the fitted priors unchanged.

### Remaining development gates for Bayesian scales (original lines 794-797)

Weak landscape effects may provide little information about scale. Partial
pooling and priors can stabilize estimation, but narrow intervals can reflect
the prior. General scale-recovery, calibration, and runtime studies remain
development gates.

## From `vignettes/bayesian-prediction.Rmd` (as of commit bfdc133e3)

### Diagnostics of the stored forecasting fit and its rejected repeat (original lines 231-241)

**How to read the fit information:** The recorded diagnostics describe the
original full chains before storing a compact set of complete joint parameter
vectors. The designated teaching fit, seed 630127, passed the unchanged raw
and derived-scale gates. The independent repeat, seed 630128, had two
divergences and remains rejected even though its Rhat and effective sample
sizes looked satisfactory. No sampler settings or thresholds were changed to
remove that rejection. Inspect every diagnostic, not Rhat alone.
This is a worked dataset, not a calibration study. The compact fit has
no full sampler for trace plots or a new chain-level convergence audit.
For your own fit, inspect `summary(forecast_fit)$diagnostics` and
`diagnose_bayes_scales(forecast_fit)` before interpreting predictions.

## From `vignettes/bayesian-computation.Rmd` (as of commit bfdc133e3)

### Evidence beyond one example: archived repeated-site simulation inputs (original lines 342-357)

On all 40 archived repeated-site simulation inputs, the nested response
proposal passed the unchanged importance checks for 19 of 20 Gaussian and
20 of 20 NB2 datasets. The conditional-coefficient alternative passed for
10 of 20 and 13 of 20, respectively. Original NUTS passed for 20 of 20 and
18 of 20. The remaining nested rejection had minimum indicator ESS 389,
below the required 400, and remains a rejection.

For the 19 Gaussian and 18 NB2 datasets where both nested correction and
NUTS passed, annual log-scale RMSE ratios were 0.999 and 1.002, and mean
95% interval-width ratios were 0.995 and 0.998. These reused-input comparisons
support this proposal choice for these tested models. They do not establish
equivalence, independent calibration, or reliability for crossed response
effects. Different eligibility sets must not be compared as if they contained
the same datasets.

### Scope of the stored NUTS eligibility flag (original lines 401-402)

The current NUTS eligibility flag checks unconstrained parameters and sampler
behavior; it does not certify annual-median Monte Carlo precision.

## From `vignettes/landscape_metric_covariates.Rmd` (as of commit bfdc133e3)

### Hard-radius metrics with the exponential-power kernel (original lines 474-479)

The native-refit backend can combine optimized landscape radii with
`kernel = "expow"`. A hard radius ignores kernel shape. Its shape and shape SE
are reported as `NA`, and its ignored shape position does not count toward K.
The likelihood remains stepped in the radius, so this support does not create
a smooth-curvature radius SE. Both FFT and direct focal raster calculation
support the categorical metric families under the same support rules.

## From `vignettes/surface_metric_covariates.Rmd` (as of commit bfdc133e3)

### Optimized hard-radius surface metrics and kernel contracts (original lines 427-438)

Optimizing a surface texture radius is well defined, but it costs more than
optimizing a kernel-weighted mean: at each candidate radius the metric is
recomputed for every sample point, since (unlike a kernel mean) it cannot be
read off a precomputed distance binning. Hard-radius metrics change when cells
enter or leave the buffer; their objectives need not be smooth. Start with fixed radii to confirm
the workflow, then optimize once the model structure is settled. Optimized
unweighted surface metrics can be combined with `kernel = "expow"` in the
native-refit backend. Their hard radius has no shape parameter, so the reported
shape and shape SE are `NA` and no shape is added to the model-selection count.
The weighted form below uses Gaussian, exponential, or exponential-power
weights. Analytic and Bayesian weighted-SQ adapters retain their narrower
Gaussian/exponential kernel contracts.

### Weighted roughness and interval coverage (original lines 450-458)

This changes both the question and the statistics. A weighted `sq` answers "at
what spatial scale does the *variability* of this surface, measured with a
smoothly decaying neighborhood, best predict the response?" Because the kernel
weights change smoothly as sigma changes (cells do not jump in and out of a hard
edge), the predictor can support analytic scale derivatives. Scale standard
errors still require eligible curvature and the selected backend's uncertainty
checks. Smoothness alone does not establish interval coverage. This is
the natural way to ask about the scale of effect of heterogeneity, as opposed to
the scale of effect of the mean.

### Connect weighted SQ to an analytic backend (original lines 502-523)

With Gaussian or exponential kernels, optimized weighted `sq` can use
`profile_glm`, `joint_nb2`, `conditional_logit`, or `joint_mixed`, subject to each
backend's model restrictions. A model can contain SQ alone or combine it with
kernel means. Keep the raw cell buffers returned by `kernel_prep()` for SQ.
Include `bin = TRUE` when the model also contains kernel means. Surface-only
inputs can omit bins, but need preparation row names that match the model rows.
Estimating other surface metrics' radii and hard-radius optimization use
`backend = "refit"`. Fixed-radius surface and landscape metrics can enter
analytic and Bayesian designs as ordinary numeric predictors. Keep the original
prepared columns and their normalization; they add no estimated bandwidth.

Live prediction and `prepare_prediction()` use the fitted normalization and
verify source specifications and row alignment. For a single-predictor NB2 GLM,
`scale_inference()` can profile the SQ scale using the existing constrained
refits. This extends the computation, not the evidence for interval coverage.
Conditional-logistic and mixed-model prediction standard errors require
`uncertainty = "conditional"`; their joint scale uncertainty remains unavailable.
The Bayesian workflow also accepts weighted SQ. It holds pooled reference
normalization fixed across years, draws and prediction rows, preserving the
meaning of coefficient priors. Its posterior intervals and diagnostics remain
different from likelihood intervals and conditional prediction SEs. Keep raw
buffers for prediction; see `vignette("bayesian-scales", package = "multiScaleR")`.

## From `vignettes/multiScaleR_Guide.Rmd` (as of commit bfdc133e3)

### Cross-validation scoring adapters still to be qualified (original lines 668-680)

See `vignette("cross-validation", package = "multiScaleR")` for worked spatial,
grouped and model-selection examples. The current scoring adapters cover
Gaussian identity, Poisson/NB2 log, binomial logit and Gamma log GLMs. NB2 theta
and Gamma dispersion are trained within each fold and frozen for assessment;
native coefficient SE conventions and nuisance/scale penalties remain separate.
Gamma and fixed-theta NB2 use native refitting. `cross_validate_clogit()`
adds the separate conditional-choice target for complete one-used-step sets,
with optional training-only spatial refitting. Its negative log scores average
over choice sets, and training selection retains native event counts. Keep
every set together in each split; use whole animals when that is the transfer
question. Generic censored survival and other model classes need
separately qualified scoring adapters. Priors, numerical guards, uncertainty
distinctions and native diagnostic requirements remain unchanged.

### Occupancy models and sample size for scale recovery (original lines 960 and 1007)

Now we'll simulate occurrence data suitable for analysis with `unmarked`. Preliminary experience has shown that simulation parameters can be harder to recover when using an occupancy model in `unmarked` with `multiScaleR`. We'll use the same raster surfaces. This section is not evaluated during vignette builds because the repeated `unmarked` refits are computationally expensive; copy the code into an interactive R session to run it.

We did a decent job of recovering simulated values. Note that the sample size was increased for this simulation. Previous explorations suggested that smaller sample sizes (e.g., 100) did not contain enough information to optimize the scales of effect (`sigma`). There is opportunity for further exploration of data needs / limitations for estimating scales of effect, especially related to more complex models such those fit with `unmarked`.

### Exponential power kernel in model selection (original line 583)

The questionably-fitting exponential power model is best-supported by AICc ranking while the Gaussian model is best-supported by BIC. Because of the added complexity (extra parameter) with the exponential power kernel, uncertain optimization, and general equivalency to other models, we should probably lean toward the simpler Gaussian kernel model. This example highlights the challenges that may be encountered when trying to parse different kernels for estimating scales of effect. From observations of performance with simulated data, use of different kernels tends to result in a similar final model and interpretation of variable effects. This, in part, is why the Gaussian kernel is the default, and likely will meet most researcher needs.

## From `vignettes/spatial-response-fields.Rmd` (source text: `vignettes/inference-contracts.Rmd` as of commit bfdc133e3)

The new vignette was written from original lines 556-643 and 841-858 of
`inference-contracts.Rmd`. The passages below are the development-status and
adapter-detail text that the vignette does not carry. The two `eval=FALSE`
code chunks (original lines 570-581 and 604-616) are not repeated here.

### What the response-field checks cover (original lines 558-568)

A native glmmTMB response field adds correlated variation to the response's
linear predictor. It does not make the ecological kernel scale vary by location.
Native refitting can retain this covariance while kernel predictors change.
Coordinate-aware prediction checks cover one shared exponential, Gaussian or
Matern field with Gaussian, log-Poisson, logit-binomial, log-NB2 or log-Gamma
responses. They use selected fixed ecological scales without a year term.
The expanded constructions fix correlation range and Matern shape while
estimating the field SD; the original Gaussian/log-Poisson exponential
references also estimate range. These checks verify native adapter behavior,
not recovery, interval calibration or general spatial model selection performance. The `joint_mixed` backend retains its
independent-random-intercept guard and does not admit these fields.

### Prediction contract and parameter counting (original lines 583-590)

The contract names the response covariance, coordinate inputs, field replication
and covariance-coordinate count. That diagnostic count describes the native
layout before parameter maps, including any coordinates fixed by the user.
Units are not guessed from factor labels. glmmTMB's likelihood degrees of
freedom supply the actual estimated-parameter penalty; mapped fixed ranges
and shapes add no parameter. Realized field values are not counted separately.
multiScaleR also counts its fitted ecological kernel parameters. A field's
correlation range and an ecological scale have different interpretations.

### Mean-prediction SEs, zero-field reference, and full-field rasters (original lines 592-602)

Native point predictions and prepared fixed-row snapshots preserve their
mean-prediction SEs, conditional on the selected ecological scales. Binomial
output is a per-trial probability even for grouped success/failure fitting
data. Poisson/NB2 maps retain native expected counts and Gaussian/Gamma maps
retain native response means. These mean SEs do not add observation variance. Setting
`re.form = NA` supplies a zero-deviation reference rather than integrating over
the field. Observation noise and ecological-scale uncertainty remain separate.
For predictions at new field locations, construct new coordinate factors with
the same projected units and prepare matching kernel inputs at those locations.
A raster's constant site scenario does not supply a spatial field at every cell.
For a full-field raster, verify the original fitting coordinates explicitly:

### Coordinate checks, batching, memory budget, and remaining gates (original lines 618-643)

New positions require the explicit `allow.new.levels = TRUE` control shown
above. The adapter fails clearly before a new-position batch enters native
prediction if that control is missing or false. Other native warnings remain
visible. Each cell supplies its own projected position rather than repeating the selected
position scenario. The shared group must match the original field. Coordinate
order, factor encoding and units are checked against retained native rows. The declared CRS must describe those original coordinates. No coordinates are
silently rescaled or transformed. Only two-dimensional projected horizontal
CRSs are accepted, including when reloading a descriptor. Geocentric,
vertical and compound definitions must be replaced explicitly with the
intended map CRS before fitting. Fit and requested SE bands remain present
when some or all cells have missing predictors. A different CRS, edited fitting context or
changed descriptor fails before prediction. Recreate the descriptor after an
intentional change.

Native calls use serial sub-batches. A conservative dense-array estimate must
fit the declared MiB budget before field expansion. This guards a planned working
set, not an operating-system memory cap or an exact accounting of native tapes.
Smaller batches reduce expansion but repeat preparation. One shared field
uses the native covariance and prediction method for the qualified families
and covariance structures above. Correlated joint analytic scale fitting,
cross-row raster covariance and spatial mixed-model CV remain separate gates.
The checked provider layout is glmmTMB 1.1.14. Mixed-model effect curves attach
the same explicit zero-field contract; their reference is not a random-field
average. GDR's native iid abundance groups remain independent effects. A
continuous GDR response field has separate normalized-target and prediction
prerequisites and is not enabled by this glmmTMB coordinate descriptor.

### Gaussian response-field cross-validation (original lines 843-858)

`cross_validate_field()` estimates a native Gaussian response field and its
correlation range using analysis rows only. Optional ecological scales and
learned transformations are also refitted in training. The held-out target is
a joint Gaussian response block at fitted parameters, with conditional field
covariance and observation noise. Coefficient, range, SD and ecological-scale
parameter uncertainty is not integrated. Native mean-prediction SEs are not
substitutes for that observation predictive covariance.

Exponential and Gaussian covariance, one shared field, constant observation SD
and ML are qualified separately from existing GLM CV. Whole field locations
remain in one split part. Buffered-out rows enter neither the model nor its
score. Compact training likelihood/count records support AIC/AICc/BIC selection;
integrated states add no penalty and an estimated range adds its native
parameter. Every failed fold/start and uncertainty distinction is retained.
Paired GDR and other non-Gaussian joint scores remain separately guarded. See
`vignette("cross-validation", package = "multiScaleR")` for the workflow.
