# Assessment of the analytic scale backends

This study checks what the public analytic estimators contribute and where
their likelihood and uncertainty results remain unresolved. Every simulated
dataset has a fixed seed and a directly calculated Gaussian kernel predictor.
The analysis stores each fit and independent reference outcome, including
optimizer failures, warnings, and unavailable uncertainty.

## Predeclared panels

| Panel | Model and design | Replication | Independent reference |
|:--|:--|:--|:--|
| `glm` | Gaussian, Poisson, grouped binomial, or NB2 with an offset | Eight datasets per family and scenario | Base R GLM likelihoods, or separate base-R NB2 likelihood optimization with the exact Poisson limit |
| `mixed` | Poisson, NB2, or zero-inflated Poisson with one site intercept and an offset | Six datasets per family and scenario | Native glmmTMB fit at each of 19 independently computed kernel values, and again at the joint fit's estimated scale |
| `bootstrap` | One-scale NB2 model evaluated at the generating scale | Four datasets per scenario, with 19 bootstrap samples per test | Profile likelihood and explicit bootstrap records from the public routine |

Each panel includes an identifiable interior scale, a weak landscape
coefficient, and a generating scale near the upper search bound. The last
scenario uses the Poisson limit for NB2 data generation. The nominal 90% Wald
interval is calculated from the reported local curvature only when a scale SE
is available. Report both the availability of an interval and its recovery of
the simulated truth. Do not treat a bound-limited fit or missing interval as
successful localization.

The profile reference checks the package's normalized likelihood at the
selected scale and whether an independently evaluated grid point has a better
likelihood. A maximum absolute same-scale difference of 0.0001 negative
log-likelihood units is a useful numerical check for an eligible reference.
A coarse grid cannot prove that an optimizer found the global maximum.

The bootstrap panel checks end-to-end reliability and Monte Carlo ambiguity
for a pointwise test. With 19 bootstrap samples, p-values have resolution
0.05. Four datasets per scenario give little information about nominal
coverage. They cannot establish bootstrap calibration or justify replacing the
package's numerical safeguards. A future calibration panel needs more
outer replicates and a substantially larger within-dataset bootstrap budget.

## Run

Run from the package source with its compiled native code available. For the
mixed panel, use ABI-compatible glmmTMB and TMB installations plus RTMB.
Use a separate output directory for each phase. The scripts resume completed
datasets and require a matching source and design fingerprint.

```text
Rscript tools/benchmarks/inference_qualification_20261112/run.R glm OUTPUT_GLM
Rscript tools/benchmarks/inference_qualification_20261112/run.R mixed OUTPUT_MIXED
Rscript tools/benchmarks/inference_qualification_20261112/run.R bootstrap OUTPUT_BOOTSTRAP
```
