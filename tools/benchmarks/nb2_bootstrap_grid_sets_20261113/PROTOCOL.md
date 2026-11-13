# NB2 bootstrap confidence-set study protocol

**Frozen before execution on November 13, 2026.** Package source starts at
`8f7a26a8ead71c29c89e725e40e61829fa11f0e7`. The earlier 24-dataset, 49-draw
pointwise study is separate and will not be pooled with these results.

## Fixed design

Use the previous study's Gaussian-kernel generator with 110 independent
rows, three fixed scenarios, and new seeds. Interior has a true scale of
240 m, kernel coefficient 0.95, and NB2 size 2. Weak uses 240 m, coefficient
0.12, and size 2. Near bound uses 730 m, coefficient 0.95, and the exact
Poisson limit. The model has one optimized bandwidth, a nuisance covariate,
and an offset. Search 50 to 800 m using `joint_nb2`, four original scale
starts (75, 200, 430, and 710 m), and `n_starts = 2` in `scale_inference()`.
Use one R process per job, no optimizer workers, and one BLAS thread.

- **Grid panel.** Simulate 16 independent datasets per scenario (48 total).
  Call public `scale_inference()` at level 0.90 with 99 bootstrap pairs at
  every point, omitting `bootstrap_values` so the package returns its full
  grid set. Supply the same fixed nodes in every scenario: 90, 140, 200,
  280, 400, 550, and 680 m. The package adds 50 m, 800 m, and the refined
  optimum. The generating scale is off the fixed grid in all scenarios.
- **Point panel.** Simulate 32 additional independent datasets per scenario
  (96 total). Test the generating scale at level 0.90 with 99 bootstrap
  pairs, using the public pointwise route. These are separate datasets and
  address pointwise inclusion, not returned set coverage.
- **Paired off-grid diagnostic.** After every full-grid call, independently
  calibrate the generating scale on the same data, using 99 pairs and its
  own prespecified seed. This checks agreement between the returned
  grid representation and a direct test of the off-grid truth. The two
  targets remain distinguishable even if their empirical rates agree.

The full-grid job runs first on its assigned process. The second process
runs the independent point panel. No case is selected for the study using
its fit or interval. All seed assignments appear in the checkpointed plan.
If the deadline approaches, both jobs finish the current dataset, retain
the pending plan, and stop. Summaries report the complete balanced prefix
and any additional attempted cases separately.

## Estimands and verification

A returned component [lower, upper] is interpreted literally as a closed
interval with endpoints equal to accepted grid values. The primary set
coverage event is membership of the generating scale in the union of those
components. A returned set that spans the search domain covers a true
scale in the domain but may provide little localization. An unavailable
set counts as no returned covering set in the unconditional operational
rate. The conditional rate uses datasets with an available set. Both
rates receive 95% exact binomial intervals. If some sets are unavailable,
also report the empirical missing-case sensitivity range obtained by
assigning all unavailable cases zero or one for coverage.

Compute set availability, empty/singleton/disconnected/full-range cases,
bound contact, total accepted length as a fraction of log search length,
largest width of any endpoint bracket on the log scale, and the interval
union's coverage. Report Jaccard overlap by log length versus the
root-refined chi-square likelihood-profile set as a descriptive comparator.
The chi-square set is not a calibration reference. Report exact
intervals for major binary outcomes and retain all per-dataset metrics.

For the coarse projection of a fine run, retain fixed nodes 50, 140, 280,
550, and 800 m and its same fitted optimum. Use the fine run's decisions
at those shared points to reconstruct the coarse components. This
comparison isolates removal of nodes from the representation. It does
not claim that a fresh coarse run with the public sequential RNG would
have identical bootstrap draws.

For the Monte Carlo comparison, reconstruct both grids using the first
49 draws at each node and compare them with the 99-draw sets. Verify
every public 99-draw decision and adjusted exact Monte Carlo interval by
independent arithmetic. Retain nominal decisions at ambiguous nodes.
At a dataset level, compute confirmed and possible membership by
reconstructing sets from the clear-accept decisions alone and from all
decisions except clear rejections. Treat any failed pair as unavailable
for that grid. This is an observed sensitivity range under per-node
ambiguity, not a simultaneous confidence band. Retain grid rows,
individual bootstrap pair records, warnings, profile audits, and failures.

With eight or more fixed points, interpolation/representation and bootstrap
error can both affect the reported set. The primary contrast is fine
versus nested coarse using identical node decisions, followed by 49 versus
99 draws on each fixed node set. The direct truth test is an additional
comparison and has its own Monte Carlo uncertainty.
