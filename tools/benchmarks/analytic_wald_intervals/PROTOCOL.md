# Analytic scale interval assessment protocol

This study compares the public analytic log-scale Wald interval with a
separately reconstructed nuisance-profile likelihood set on the same datasets.
The one Gaussian kernel predictor is standardized at each candidate scale. A
fixed nuisance covariate and known offset are present in every data-generating
and fitting model. Native analytic defaults determine optimizer starts, Wald
availability and interval construction.

## Design

Four independent-response families are simulated. Gaussian observations
have residual SD 0.65. Poisson observations have a log mean. Grouped binomial
observations have five trials per row and a logit probability. NB2 observations
have log mean and size two, except in the separately simulated exact Poisson
dispersion face. Both analytic NB2 backends fit the same NB2 datasets.

| Scenario | Generating scale | Landscape coefficient | Observation variance |
| --- | ---: | ---: | --- |
| Interior | 240 | 0.95 | Stated family |
| Weak | 240 | 0.12 | Stated family |
| Near upper bound | 730 | 0.95 | Stated family |
| Exact upper bound | 800 | 0.95 | Stated family |
| Exact lower bound | 50 | 0.95 | Stated family |
| NB2 Poisson face | 240 | 0.95 | Poisson limit |

Each family-scenario combination has 48 fixed-seed datasets at each of
(n = 75) and (n = 240). The scale search spans 50 to 800 map units.
The cohort has 2,016 datasets and 2,592 backend fits. The exact Poisson
dispersion face contributes NB2 datasets only. The near and exact scale
boundaries retain finite NB2 dispersion to distinguish those regimes.

## Matched inference and numerical checks

The reference recomputes the Gaussian kernel directly from the distance bins
at each scale. It refits Gaussian, Poisson and binomial coefficients through
base R GLMs and calculates normalized response densities explicitly. NB2
uses direct R densities, Fisher scoring for coefficients, and dispersion
profiling over the same admissible range as its paired backend. The
`joint_nb2` reference includes alpha = 0 exactly. The `profile_glm` reference
uses the same finite theta ceiling as the package.

The unrestricted reference refines every detected local minimum and both
end intervals. A second, denser 81-point grid checks the original 41-point
grid and both include the true and package-estimated scales. Profile set
endpoints are root refined on the log scale. Each distinct component is
retained along with its boundary flags. Grid convergence is a resolution
check. It cannot establish that narrower components are absent. The nominal
reference is the chi-square likelihood-ratio set, at both 90% and 95%.
The study uses its inclusion of the known generating scale as the primary
coverage decision and checks that the geometric union gives the same answer.

The public default `summary()` Wald interval is checked against the package
covariance formula at 95%. The 90% calculation uses the same estimate,
uncertainty availability and log-scale covariance with the 90% quantile.
For the first replicate of each design cell, also compare the public 95%
profile interval with the corresponding component of the independent set.

Preserve all planned datasets, fit messages, analytic starts, numerical
reference errors, warning records, and detected profile components. Report
coverage among returned intervals, unconditional availability and coverage,
matched comparisons, boundary/full-range frequencies, and Wilson 95%
Monte Carlo intervals. A missing interval contributes no success to the
unconditional measure. A bound-limited or full-range set may include the
truth and may provide limited localization. At exact scale or dispersion
bounds, the ordinary chi-square cutoff remains a procedure to assess.

## Numerical follow-ups

The boundary-start check applies to default backend fits with an independent
global likelihood improvement above 1e-5 NLL units. It uses eight fixed starts
including both scale bounds. The follow-up records original and new likelihoods,
the independently reconstructed optimum, and all starting points. It leaves the
primary estimates and denominators as observed.

Early coverage results motivated a post hoc curvature check. Before running
the check, we fixed the subset to the first four replicates of each cell.
We independently differentiate the nuisance-profile NLL when the default Wald
interval and the reference profile set are available. Use a five-point central stencil on the log scale
and repeat with half the step. Report numerical availability, relative change
after step halving, and relative error between the inverse-curvature SE and
the package's log-scale SE. Keep the locally estimated scale in both
calculations. Local curvature agreement and empirical coverage are separate
comparisons.

## Reproduce

Run the two shards from the package source root with the matching compiled
package DLL available. Set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and
`MKL_NUM_THREADS` to 1 before starting. Use a new output directory for a new design or source
fingerprint. The run refuses to resume with different inputs.

```text
Rscript --vanilla tools/benchmarks/analytic_wald_intervals/run.R 1 RUN_DIRECTORY
Rscript --vanilla tools/benchmarks/analytic_wald_intervals/run.R 2 RUN_DIRECTORY
Rscript --vanilla tools/benchmarks/analytic_wald_intervals/summarize.R RUN_DIRECTORY OUTPUT_DIRECTORY
Rscript --vanilla tools/benchmarks/analytic_wald_intervals/check_boundary_starts.R OUTPUT_DIRECTORY/metrics.csv OUTPUT_DIRECTORY
Rscript --vanilla tools/benchmarks/analytic_wald_intervals/check_local_curvature.R OUTPUT_DIRECTORY/metrics.csv OUTPUT_DIRECTORY
```

The two shards use disjoint seeds within every cell. The dataset records in
the run directory are atomic and include all detected optimization starts
and independent reference warnings or failures. The generated CSVs are
derived from these records.
