# Analytic Wald scale intervals and full-range likelihood sets

The matched study completed 2,016 datasets and 2,592 backend attempts. It
compared the public log-scale Wald interval with independently reconstructed
nuisance-profile likelihood sets on the same datasets. The two NB2 backends
shared every NB2 dataset and differed only in their allowed dispersion range.
The design used one Gaussian kernel mean, a fixed covariate, a known offset,
and 48 datasets per response, scenario, and sample size. Sample sizes were 75
and 240. See the [protocol](PROTOCOL.md) for the full design and commands.

## Coverage and availability

The table reports nominal 95% results by scenario. Each cell gives counts
over backend attempts. NB2 contributes two paired rows per dataset to these
scenario totals. The denominators in the two conditional coverage columns
include only results that returned the specified interval or set. Coverage
over all attempts counts each missing interval or set as unavailable and
uncovered.

| Scenario | Attempts | Wald returned | Wald covered among returns | Wald covered over all | Profile returned | Profile covered among returns | Profile covered over all |
|:---|---:|---:|---:|---:|---:|---:|---:|
| Interior | 480 | 476 | 401/476 | 401/480 | 474 | 442/474 | 442/480 |
| Weak effect | 480 | 276 | 243/276 | 243/480 | 480 | 452/480 | 452/480 |
| Near upper bound | 480 | 323 | 280/323 | 280/480 | 470 | 453/470 | 453/480 |
| Exact upper bound | 480 | 240 | 199/240 | 199/480 | 472 | 450/472 | 450/480 |
| Exact lower bound | 480 | 252 | 211/252 | 211/480 | 480 | 458/480 | 458/480 |
| NB2 Poisson dispersion face | 192 | 0 | unavailable | 0/192 | 192 | 184/192 | 184/192 |

The nominal 95% interior results by family are shown below. Brackets give
Wilson 95% Monte Carlo intervals for coverage over all attempted datasets.
The two NB2 backends had the same coverage counts. In the 470 interior
backend attempts with both results available, both covered in 395 cases,
only the profile covered in 43, and neither covered in 32.

| Response | Wald available | Wald covered over all [95% interval] | Profile available | Profile covered over all [95% interval] |
|:---|---:|---:|---:|---:|
| Gaussian | 96/96 | 88/96 [0.844, 0.957] | 96/96 | 92/96 [0.898, 0.984] |
| Poisson | 96/96 | 82/96 [0.770, 0.911] | 96/96 | 88/96 [0.844, 0.957] |
| Grouped binomial | 96/96 | 85/96 [0.806, 0.935] | 96/96 | 90/96 [0.870, 0.971] |
| NB2, each backend | 94/96 | 73/96 [0.666, 0.835] | 93/96 | 86/96 [0.819, 0.942] |

Interior NB2 covered in 35/48 and 38/48 Wald attempts at the two sample sizes.
The corresponding profile sets covered 42/48 and 44/48. At nominal 90%,
interior Wald covered 371/480 backend attempts and profiles covered 406/480.
[Cell results](results/cell-rates.csv) retain both levels and all denominators.
[Matched NB2 records](results/nb2-matched-backends.csv) show the paired
decisions. A single weak NB2 dataset at (n=75) returned a Wald interval
under `profile_glm` only.

## Localization and numerical checks

At 95%, 258 of 480 weak-effect sets covered the entire search range, and eight
had disconnected regions. Profiles reached a bound in 425 of 470 near-upper
returns. Across the cohort, 746 of 1,567 returned Wald intervals crossed a
search bound. Among the 740 of those attempts with an available profile,
739 profile sets reached a bound. A further 148 profiles reached a bound
even though both Wald endpoints stayed inside. Neither endpoint pattern
establishes coverage or proves that a scale is localized.

Both independent profile grids gave the same number of detected components
for every available set. The largest absolute log-endpoint difference after
doubling grid resolution was 1.70e-8. Membership of the generating scale
in the reconstructed union agreed with its pointwise likelihood-ratio
decision throughout. The public `summary(fit, profile = TRUE)` returned the
connected region in all 54 predeclared checks. Its largest endpoint difference
from the corresponding independent region was 0.00139 map units.

The maximum absolute same-scale negative log-likelihood difference was
3.67e-7 across available paired calculations. The only two global NLL
discrepancies above 1e-5 were the two backends on one weak NB2 dataset.
Their default starts all passed at a scale near 135.70. The corresponding
default NLL was 410.437503800088. Both fits returned Wald intervals.
Adding starts at both bounds and six fixed interior scales found a better
fit at the upper bound of 800. The new NLL was 410.430674036639, within
9.7e-13 of the independent optimum. Both refits withheld the Wald interval
for the bound and competing modes. [Boundary-start follow-up](results/boundary-start-diagnostics.csv)

Early coverage results prompted a separate curvature check. This compared
the package log-scale SE with the inverse-curvature SE from the independent
nuisance profile at the package estimate. Of 216 planned backend records,
140 had both results, including all 40 interior records. In 74 records
Wald was unavailable. Two further records had an unavailable reference.
The maximum relative SE discrepancy was 9.78e-6, including the interior
subset. Halving the finite-difference step changed curvature by at most
4.04e-8 in relative units. The observed interior coverage deficit therefore
persists in a design where the local curvature calculation agrees with an
independent calculation. [Curvature follow-up](results/local-curvature.csv)

## Failures and scope

The original records include 2,590 successful package fits and 2,568
independent profile returns. Two `joint_nb2` attempts at the Poisson face and (n=75)
failed with a nonfinite value during optimization. Twelve NB2 datasets
failed the independent reference coefficient iteration limit in both
backends. Their 24 reference failures are preserved, with no profile or
coverage result assigned. These failures occurred among the interior,
near-upper, and exact-upper scenarios. [All attempts](results/metrics.csv) ·
[Reference failures and warnings](results/reference_warnings.csv.gz)

The experiment addresses one optimized scale in independent GLMs with a
Gaussian kernel mean. Its numerical agreement applies to calculations for
which both methods returned. The coverage denominator retains every attempt. Other kernels, weighted surface metrics, multiple scales,
dependence, and calibration at statistical boundaries require their own
assessments. The chi-square profile cutoff is approximate, and a full-range
set gives little information about the location of the scale.

For a scale conclusion, inspect likelihood across the admissible range.
For supported one-scale NB2 fits, `scale_inference(fit)` retains the detected
regions and records when a region reaches a search bound. If the profile
suggests competing optima, compare starts near those regions and at the
bounds. In `profile_control$starts`, provide distances divided by
`unit_conv`. The default Wald interval and fitting rules remain unchanged.
