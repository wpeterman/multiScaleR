# Reproduce the recorded study tables

The archived runner outputs can be summarized in any directory with R,
without loading multiScaleR or refitting. From the package root, run:

```text
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/summarize.R tools/benchmarks/nb2_bootstrap_grid_sets_20261113/archive/grid tools/benchmarks/nb2_bootstrap_grid_sets_20261113/archive/point OUTPUT tools/benchmarks/nb2_bootstrap_grid_sets_20261113/archive/grid-a tools/benchmarks/nb2_bootstrap_grid_sets_20261113/archive/grid-b
```

See the [pointwise report](POINTWISE_RESULTS.md) for the separate
96-dataset panel. The grid panel retains every returned component and
each bootstrap pair from the full-grid call in the per-dataset checkpoints.
The configuration, source hashes, session records, and scheduling split
are under `archive/`. Verify the study payload against `SHA256SUMS`. For an
interrupted run, the checkpoint manifest and plan identify pending
datasets. The summary uses completed datasets in its analysis
denominators and reports the complete balanced prefix separately.

After summarizing, generate the finite-grid report from the archived
tables without rerunning either simulation.

```text
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/grid_report.R OUTPUT RESULTS.md
```
