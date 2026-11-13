# Joint radius profile numerical assessment

The benchmark evaluates point optimization at a declared hard-radius grid. It asks whether the package reaches the same likelihoods and chooses the same grid row as an independent dense native search. It does not estimate interval coverage.

## Design and reference

Three prespecified seeds are used in each of four scenarios. Every dataset has 60 observations and two hard-radius predictors, Shannon diversity and unweighted terrain standard deviation. Each fit attempts the same six radius combinations. A 40 m diversity radius contains a single cell, which makes the diversity predictor constant for two of the six rows. These rows remain in the results. The other four combinations are fitted in every replicate.

| Scenario | Smooth predictor and response |
|:--|:--|
| Gaussian interior | Gaussian kernel mean generated with a 104 m bandwidth within the 40 to 240 m fit range |
| Gaussian lower bound | Gaussian kernel mean generated with a 15 m bandwidth below the fitted range |
| Gaussian flat | One habitat value per observation, repeated across all distances to make the bandwidth likelihood flat |
| Conditional choice | One used alternative per set of three, with an exponential-power kernel mean generated at 104 m and shape 2.5 |

Reference covariates are built directly from cell values and distances. The occupied production bins each contain one distinct distance, so binning preserves the exact raw-cell weight. The Gaussian reference uses QR least squares and the maximized Gaussian likelihood from the residual sum of squares. The choice reference uses the native C Cox fit with one event per stratum. Its likelihood is checked against an exact native formula fit before every replicate. The Gaussian scan uses 121 log-spaced bandwidths, including both bounds, and refines all bracketed local maxima. The choice scan uses 41 log-spaced bandwidths by 41 log-spaced shapes from 0.75 to 50, including bounds, and refines the dense leaders and the best point from each of five shape bands. The independent scan and its refinements supply a numerical reference on the chosen bounds. They do not prove a global maximum for every landscape.

The prespecified numerical agreement threshold is an absolute log-likelihood difference of `1e-5` per scored hard-radius combination. The candidate table also evaluates the independent likelihood at the package's returned bandwidth and shape, which separates likelihood reconstruction error from a difference in the optimized point.

## Results

The package and independent search scored the same 48 of 72 hard-radius
combinations. All 24 rejected combinations had zero diversity variance.
The independent search evaluated 24,528 dense points for the other
combinations without a failed likelihood evaluation. All 144 package
smooth starts at those combinations passed the eligibility check. The
package and independent search selected the same hard-radius row in all
12 datasets.

| Scenario | Compared combinations | Largest absolute optimization gap | Selected smooth bandwidths (m) | Conditional curvature |
|:--|--:|--:|:--|:--|
| Gaussian interior | 12 | `9.33e-14` | 101.46 to 104.57 | Available in three of three |
| Gaussian lower bound | 12 | `2.13e-14` | 40 in all three | SE withheld at the bound |
| Gaussian flat | 12 | `3.29e-14` | 120 to 153.35 | SE withheld for unresolved curvature |
| Conditional choice | 12 | `2.16e-9` | 74.23 to 159.31 | Available in three of three |

The largest likelihood reconstruction difference at package-reported
parameters was `3.55e-14`. The largest difference among selected model
likelihoods was `9.33e-14`. The maximum optimization gap occurred at a
nonselected choice grid row. Its magnitude is below the prespecified
tolerance. The choice fits selected shapes between 1.27 and 5.23 across
the three datasets.

In the flat scenario, the largest range across the independent dense
bandwidth scan was `6.75e-14`. The returned bandwidths are arbitrary points
on that flat profile. Their conditional Wald SEs are unavailable. In the
lower-bound scenario the likelihood is maximized at the fitted lower
limit, so the returned 40 m is a constrained optimum. The study supports
point optimization on these declared grids. It does not assess hard-radius
confidence sets or coverage of the available conditional smooth intervals.

## Issues found and retained

The first run rejected all 36 Gaussian grid rows with valid reference
objectives. Native data recovery changed response names and removed
`scale()` provenance attributes. The fitted row IDs, responses, weights,
and offsets were equal. The original run and the intermediate names-only
attempt are retained in `before-frame-fix/` and `after-name-fix/`.
The audit now ignores names and those two provenance attributes while
checking exact row order, values and response structure.

After that correction, two flat fits returned finite conditional
bandwidth SEs of approximately 515,000 and 728,000 m. These are preserved
in `before-curvature-fix/`. The conditional Hessian now compares its
smallest curvature in log coordinates with a finite-difference resolution
floor based on machine precision, the objective magnitude, and the
smallest local step. It withholds the smooth Wald block if curvature
cannot be resolved. The likelihood table and selected point are still
available.

## Reproduce

Run from the package root in the worktree containing the joint profile implementation and a compiled package library.

```sh
Rscript tools/benchmarks/radius_profile_20261112/run.R
```

The optional first argument writes a copy of the output to another directory. The output includes every hard-radius combination in `candidates.csv`, attempted package smooth starts in `starts.csv`, selected fits and uncertainty flags in `replicates.csv`, and the R session and source hashes in `session.txt`.
