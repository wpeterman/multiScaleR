# Suggested table-reading sequence

Start with `datasets.csv` for the planned, completed, and pending sample.
The `prefix` flag identifies complete waves with equal replication across
the three scenarios. Compare the two denominator scopes before using a
rate from an interrupted run.

Then use `set_exact_binomial.csv` for the outcome counts and exact
intervals. `unconditional_coverage` counts an unavailable set as failure
to return a set covering truth. `conditional_coverage` uses available
sets as its denominator. The `zero_log_length` event counts available
sets with no positive-width component. The corresponding Monte Carlo
sensitivity bounds are in `set_summary.csv`. For an individual dataset,
see `returned_components.csv` and `grid_nodes.csv`.

Compare coarse and fine results with identical decisions at the shared
nodes before comparing 49 and 99 draws. `paired_comparisons.csv` records
each change in coverage and component endpoints, and `comparison_summary.csv`
adds exact intervals for the fraction of cases that change and records
both directions of change. It also reports median log-length Jaccard
overlap.

Use `point_exact_binomial.csv` for the independent pointwise panel and the
paired off-grid diagnostic. The `panel` column identifies the population
of datasets. The two panels address different targets and have different
seeds, so their outcomes are not pooled. The file
`point_chi_square_comparison.csv` retains the chi-squared pointwise
decision and both directions of disagreement on the same datasets.
