# joint_mixed log-scale Wald intervals for kernel sigma

This study measures how often the new `joint_mixed` 95% Wald interval for
sigma contains the generating value. The standard error comes from the TMB
`sdreport` of the selected start: the square root of the `kernel_par` entry of
`cov.fixed`, the inverse Hessian of the Laplace-approximated marginal
likelihood. The interval is `sigma x exp(+/- 1.96 x SE)` on the log scale, the
same form the analytic GLM backends report. Coefficient and prediction
uncertainty were not changed and are not assessed here.

## Design

The design was fixed before any result was seen ([run.R](run.R),
[design.R](design.R)). Each dataset has 200 points in 40 sites, five points per
site, with a normal random intercept (SD 0.5). One Gaussian kernel mean enters
the log mean with coefficient 0.8, alongside a fixed covariate (coefficient 0.3)
and a known offset. The binned radial landscape is the generator of
[analytic_wald_intervals](../analytic_wald_intervals/REPORT.md). The search
range is 50 to 800 map units.

| Factor | Levels |
|:---|:---|
| Response | Poisson; NB2 with size 2 |
| Generating sigma | 120, 240, 450 map units |
| Replicates | 150 fixed-seed datasets per response and sigma (900 in total) |
| Fitting | `multiScale_optim(backend = "joint_mixed")` with its six default starts |

For a likelihood-ratio reference, each dataset also refits the native
`glmmTMB` model with the covariate held at the generating sigma. Twice the
difference between that negative log-likelihood and the joint optimum is
compared with the chi-square cutoff. This checks whether a profile interval
would contain the truth without computing its endpoints.

## Results

All 900 fits returned, none raised a warning, and no estimate sat at a search
bound. Every reported limit matched `summary(fit)` to 1e-10. One NB2 SE was
withheld because a competing start within the 95% likelihood threshold moved
sigma by more than 0.5 log units. No log-scale SE exceeded one.

Brackets give Wilson 95% Monte Carlo intervals. "Reported" uses the shipped
rule. With no bounded fits, it differs from using every interior SE only by
the one withheld NB2 interval.

| Response | Sigma | Reported | Wald covered | Wald coverage | LR covered | SD of log error | RMS log SE |
|:---|---:|---:|---:|:---|---:|---:|---:|
| Poisson | 120 | 150 | 141 | 0.940 [0.890, 0.968] | 141/150 | 0.063 | 0.058 |
| Poisson | 240 | 150 | 132 | 0.880 [0.818, 0.923] | 140/150 | 0.123 | 0.116 |
| Poisson | 450 | 150 | 144 | 0.960 [0.915, 0.982] | 142/150 | 0.059 | 0.062 |
| NB2 | 120 | 150 | 145 | 0.967 [0.924, 0.986] | 144/150 | 0.086 | 0.092 |
| NB2 | 240 | 149 | 115 | 0.772 [0.698, 0.832] | 133/150 | 0.197 | 0.159 |
| NB2 | 450 | 150 | 144 | 0.960 [0.915, 0.982] | 142/150 | 0.102 | 0.097 |
| All | | 899 | 821 | 0.913 [0.893, 0.930] | 842/900 | 0.115 | 0.103 |

At 90%, Wald intervals covered 780 of 900 interior fits and likelihood-ratio
decisions covered 792. [Per-scenario counts](results/scenario-summary.csv)
and [per-dataset records](results/metrics.csv) hold every denominator.

## Where coverage falls short

At sigma 120 and 450, Wald coverage was close to nominal in both families
(141 to 145 of 150). The standardized errors had an SD near one, and larger
SEs went with larger errors.

At sigma 240 the interval undercovered: 132 of 150 for Poisson and 115 of 149
for NB2. The misses were not bias. Median log errors were 0.016 and -0.003,
and NB2 missed equally often above and below (17 and 17). The problem is that
the local curvature sometimes overstates precision. For NB2 at 240, the SE was
negatively correlated with the absolute log error (r = -0.49). Coverage was
0.42 in the third of fits with the smallest SEs, against 0.90 and 1.00 in the
other two thirds. Fourteen NB2 fits had errors more than three SEs from the
truth. In 17 NB2 and 8 Poisson datasets the Wald interval missed while the
likelihood-ratio decision covered. The reverse never happened at sigma 240.
This pattern fits a likelihood that is sharply peaked near the optimum but
stays within the cutoff over a wider range of scales. It depends on how fast
the kernel covariate changes with sigma in a given landscape, not on the
mixed model.

The likelihood-ratio reference also fell short at NB2 sigma 240 (133 of 150),
so part of the deficit is the small-sample behavior of the scale estimate
itself.

The pattern matches the analytic GLM study, which used the same landscape
generator at sigma 240. There, interior `profile_glm` and `joint_nb2` Wald
intervals covered 401 of 476 (0.84), with Poisson 82 of 96 and NB2 73 of 96
per backend. The `joint_mixed` SEs are no less reliable than the analytic
GLM SEs the package already reports, and they share the same weakness.

## Recommendation

Ship the SEs with a caution in the documentation, worded like the existing
`profile_glm` caution. The calculation is correct: it reproduces the
`sdreport` covariance and matches the public limits, and coverage is near
nominal in two of the three scales. But a returned interval can undercover
badly in some landscapes even when the fit is interior and every check
passes. In this study, 95% intervals held the truth in 77% to 97% of datasets,
depending on the scale and response. Users who need a firm interval should
inspect the likelihood across the search range. Unlike the GLM backends,
`joint_mixed` has no profile route yet (`profile_sigma()` and
`summary(fit, profile = TRUE)` reject it). The likelihood-ratio decisions here
(842 of 900 covered) suggest that a profile interval for `joint_mixed` would be
worth adding.

## Scope

One Gaussian kernel mean, one random intercept, 200 points, and one landscape
generator. No fit reached a search bound, so the bound rule was exercised only
by the unit tests and a bounded integration fit, not by this study. Exponential
kernels, weighted surface SQ, several scales, zero inflation, fewer groups,
and weaker effects need their own assessments. The run took 89.6 minutes of
wall time on 10 workers (730 fit-minutes in total).

To reproduce, from the source root:

```
Rscript tools/benchmarks/joint_mixed_wald/run.R <output_dir> 10
Rscript tools/benchmarks/joint_mixed_wald/summarize.R <output_dir>
```
