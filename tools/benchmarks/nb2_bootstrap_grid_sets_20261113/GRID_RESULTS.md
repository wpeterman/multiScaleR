# NB2 bootstrap finite-grid confidence-set results

The grid panel completed 48 of 48 planned datasets. The balanced prefix contains 48 datasets.
It attempted 46,332 grid refit pairs (0 failed, 13 warned) and 4,752 separately seeded truth-test pairs (0 failed, 0 warned).
The balanced prefix contains 35 ambiguous 99-draw node decisions among 468 tested nodes.

## Coverage and availability

Coverage means membership of the generating scale in the literal union of returned closed intervals.
The unconditional rate counts an unavailable set as no returned covering set.
The conditional rate uses only available sets.
Parentheses give exact 95% binomial intervals across datasets.
The direct test column gives pointwise inclusion among datasets with an available direct decision.

| Scenario | Available | Unconditional coverage | Conditional coverage | Direct truth inclusion |
|:--|--:|--:|--:|--:|
| Interior | 16/16 | 6/16 (0.152 to 0.646) | 6/16 (0.152 to 0.646) | 13/16 (0.544 to 0.960) |
| Weak effect | 16/16 | 15/16 (0.698 to 0.998) | 15/16 (0.698 to 0.998) | 15/16 (0.698 to 0.998) |
| Near upper bound, Poisson limit | 16/16 | 12/16 (0.476 to 0.927) | 12/16 (0.476 to 0.927) | 15/16 (0.698 to 0.998) |

The fine grid and direct truth test disagreed in 10 balanced-prefix datasets with both decisions available.
The generating scale was off the fixed grid. Its direct test also used a separately seeded bootstrap.
The comparison reflects both differences in the target and Monte Carlo uncertainty.

## Localization

Lengths are measured on the log-scale domain. Span is a fraction of the full search length.
Jaccard is descriptive overlap with the root-refined chi-squared profile set.
The latter is a likelihood comparator. It supplies no calibration truth.

| Scenario | Full range | Bound contact | Disconnected | Singleton component | Zero total length | Median span | Median profile overlap | Maximum log bracket width |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| Interior | 0/16 | 0/16 | 0/16 | 0/16 | 0/16 | 0.045 | 0.277 | 0.357 |
| Weak effect | 11/16 | 13/16 | 1/16 | 0/16 | 0/16 | 1.000 | 1.000 | 0.588 |
| Near upper bound, Poisson limit | 0/16 | 14/16 | 0/16 | 2/16 | 2/16 | 0.059 | 0.558 | 0.212 |

Available zero-length sets by method

- `fine99` had near_bound-02, near_bound-14.
- `coarse99` had near_bound-01, near_bound-02, near_bound-03, interior-05, near_bound-08, interior-09, weak-09, near_bound-09, near_bound-10, near_bound-11, interior-14, near_bound-14, near_bound-16.
- `fine49` had near_bound-02, near_bound-14.
- `coarse49` had interior-01, near_bound-01, near_bound-02, near_bound-03, interior-05, near_bound-08, interior-09, weak-09, near_bound-09, near_bound-10, near_bound-11, interior-14, near_bound-14, near_bound-16.

For `near_bound-02`, the accepted component was [800, 800] m. The adjacent rejected node was 680 m. This singleton contains only 800 m, so the 730 m generating scale was outside it.

Interior sets had these accepted-node counts by dataset: 2 nodes in 9 datasets, 3 nodes in 7 datasets. 16 of 16 sets had exactly one component.
Both fixed nodes flanking 240 m (200 and 280 m) were accepted in 6 covering and 0 noncovering datasets. Among 10 noncovering sets, 3 lay entirely below 240 m and 7 entirely above it.

## Grid and Monte Carlo sensitivity

The coarse projection removes nodes without refitting. Its decisions at shared nodes are identical to those of the fine run.
The 49-draw result uses the first 49 bootstrap draws at each node. It is paired with the 99-draw result.
The last column bounds observed coverage when ambiguous node decisions and unavailable sets can vary.
It is a sensitivity range for these datasets.

| Scenario | Fine 99 | Coarse 99 | Fine 49 | Coarse 49 | Fine 99 cases with ambiguity | Fine 99 sensitivity range |
|:--|--:|--:|--:|--:|--:|--:|
| Interior | 6/16 | 1/16 | 5/16 | 0/16 | 5 | 0.250 to 0.562 |
| Weak effect | 15/16 | 14/16 | 14/16 | 14/16 | 6 | 0.875 to 1.000 |
| Near upper bound, Poisson limit | 12/16 | 5/16 | 12/16 | 5/16 | 7 | 0.500 to 0.812 |

The paired table reports the fraction of datasets whose coverage event changes. Parentheses give exact 95% intervals.
The final two counts in each cell give datasets covered only by the first and only by the second method, respectively.
For the direct test, the paired event is truth inclusion.

| Comparison | Interior | Weak effect | Near bound |
|:--|--:|--:|--:|
| coarse99_vs_fine99 | 5/16 (0.110 to 0.587), 0 / 5 | 1/16 (0.002 to 0.302), 0 / 1 | 7/16 (0.198 to 0.701), 0 / 7 |
| fine49_vs_fine99 | 1/16 (0.002 to 0.302), 0 / 1 | 1/16 (0.002 to 0.302), 0 / 1 | 0/16 (0.000 to 0.206), 0 / 0 |
| coarse49_vs_coarse99 | 1/16 (0.002 to 0.302), 0 / 1 | 0/16 (0.000 to 0.206), 0 / 0 | 0/16 (0.000 to 0.206), 0 / 0 |
| fine99_vs_direct99 | 7/16 (0.198 to 0.701), 0 / 7 | 0/16 (0.000 to 0.206), 0 / 0 | 3/16 (0.040 to 0.456), 0 / 3 |

Calibration remains experimental. The evidence is specific to the declared grids, a single Gaussian kernel predictor, and a finite bootstrap budget.
See the [separate pointwise panel](POINTWISE_RESULTS.md) for the 96-dataset truth-test results.
