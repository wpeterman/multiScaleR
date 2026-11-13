# NB2 bootstrap confidence-set assessment

The recorded study has a [finite-grid report](GRID_RESULTS.md) and an
[independent pointwise report](POINTWISE_RESULTS.md). [Archive instructions](ARCHIVE.md)
show how to regenerate the tables without refitting.

See [the frozen protocol](PROTOCOL.md) for the sampling design, analysis
targets, and run limits. The integrated worktree must have its compiled
native library available. The study calls the public package fit and
inference functions with their existing convergence checks.

```text
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/run.R grid GRID_OUTPUT
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/run.R point POINT_OUTPUT
```

Run at most one process for each panel. For each process, set
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `MKL_NUM_THREADS` to 1.
A completed dataset writes an atomic, compact RDS checkpoint containing
the observed grid, bootstrap pair records and optimizer audits. A resumed
run verifies the saved plan and source hashes before fitting pending cases.
An empty `STOP` file in the relevant output directory requests a stop
after the current dataset. Remove it before resuming.

Summarize both outputs without refitting.

```text
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/summarize.R GRID_OUTPUT POINT_OUTPUT REPORT_OUTPUT
```

`datasets.csv` contains the complete plan, checkpoint status and refit
counts. `sets.csv`, `points.csv`, and `paired_comparisons.csv` retain
dataset-level analysis results. The exact binomial tables give the
event count and its denominator for each rate. Every summary includes
the complete balanced prefix and the full set of attempted cases as
separate scopes. See [diagnostic interpretation](READOUT.md).

If the point panel finishes while the grid panel is still running, stop
the original grid runner after its current dataset. The remaining grid
IDs can then be allocated to two disjoint output directories without
changing the original plan, inputs, or seeds.

```text
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/allocate.R GRID_OUTPUT ALLOCATION
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/shard.R GRID_OUTPUT GRID_A ALLOCATION/a.txt
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/shard.R GRID_OUTPUT GRID_B ALLOCATION/b.txt
Rscript tools/benchmarks/nb2_bootstrap_grid_sets_20261113/summarize.R GRID_OUTPUT POINT_OUTPUT REPORT_OUTPUT GRID_A GRID_B
```

The shard driver verifies the original source and runner hashes, and the
summary rejects overlapping assignments. Maintain at most two
single-thread processes overall.
