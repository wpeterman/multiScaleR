# Native GDR spatial response contract, 0.7.31

`fit_gdr_response_field()` is an opt-in compiled native TMB route returning an
identified `unmarkedFitGDRField`. The initial scope is Poisson abundance,
half-normal detection, density per square kilometer, one primary period and
complete canonical 75-meter point-circle bins. Supply a prespecified exponential
range and one consistent projected location per original abundance state.
One Gaussian abundance SD is estimated; the optional original independent
distance-detection observer SD is retained. No year field or population prior
is introduced. The private RTMB fitting admission remains false.

The native finite-K distance/removal likelihood and normalized Gaussian SDs
are preserved. Correlation-whitened abundance coordinates use the transformed
native random design. Both original and independent cold coefficient/SD starts
use zero latent starts and must pass the unchanged convergence, score,
current-mode latent precision, full fixed/joint curvature and agreement gates.
Every attempt and failure has a compact write-once audit. No model, tape,
chain, covariance array or candidate array is checkpointed in those audits.

Native coefficient and SD counts remain. The fixed range and integrated
realized field states add zero penalty parameters. Estimated ecological scales
and shapes keep their ordinary multiScaleR counts. The six-parameter observer
fixture and five-parameter no-observer fixture verify this distinction.

`make_gdr_refitter()` retains frozen component design recipes and the original
field/protocol through raw-reference scale refits. Scale fitting is serial:
refit backend, one core and PSOCK=FALSE. Unqualified native fit/tape transport
is guarded. Prediction/refit closure environments retain only their declared
numeric recipes, protocol data and namespace parent, without an unused original
fit. Native caller thread settings are restored after each attempt.

Point and snapshot predictions retain native component labels. Conditional
density predictions at new projected coordinates combine full fitted
coefficient/field cross-covariance with separately identified new-field
residual variance. Neither is observation noise or ecological-scale uncertainty.
Inverse-link point predictions are not lognormal field averages. re.form=NA
sets response effects to zero; it does not integrate them. Native covariance
and normal intervals use retained physical field coefficients/covariance and
do not rerun a native TMB tape.

Actual Gaussian kernel refitting, numeric snapshots, zero-effect curves and
bounded coordinate density/SE maps are tested. Logical raster columns retain
their declared 0/1 types and invalid-encoding guards. Maps require matching planar
CRS, one worker and bounded batches. Pure random-group labels are not automatic
fixed raster scenarios. Factors that also enter fixed components retain their
explicit fitted-level scenario requirement. Existing model routes are unchanged.
No full-grid covariance is retained. Numeric fresh-worker predictions use only
the portable field packet; native fitting workers remain guarded.

Canonical code/formals and captured numeric projection inputs are bound
independently of runtime compilation. Actual code, controls, native layouts,
row identities, geometry, fitted target, covariance and diagnostic mutations
are rejected. Provider/version/DLL and source/archive identities remain explicit.

Estimated ranges, broader GDR mixtures/keys, annual availability, spatial GDR
CV, field simulation, latent-abundance posterior summaries and profile intervals
require separate qualification. Private optimized RTMB/native agreement checks
are derivative/implementation references, not public RTMB fitting admission,
an assessment of Laplace approximation error, empirical recovery or calibration.
Spatially varying ecological kernels remain a separate guarded development.

The existing S3 covariance method has an explicit stats::vcov registration
beside the new S4 methods. A fresh worker requires the package namespace
without attaching it and verifies registration and native conditional
covariance. Ordinary covariance semantics and every method body are unchanged.
