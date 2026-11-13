# Replicated NB2 pointwise bootstrap assessment

## Target

The study estimates how often the package's finite-bootstrap test includes the generating scale at a nominal 90% confidence level. For each dataset, `scale_inference()` calibrates the likelihood-ratio test at the generating scale alone. The likelihood is maximized under the generating scale and unrestricted scale, with coefficients and NB2 dispersion reestimated in both fits. Each bootstrap draw has a matched constrained and unrestricted refit. The generating law for the draw is the fitted constrained model.

Pointwise inclusion is the event needed to invert a continuum of level tests. This run evaluates only one scale per dataset. It cannot determine coverage, localization, or topology of the package's returned finite-grid confidence sets. Those also depend on the grid and the other pointwise tests.

## Design

Each scenario has 110 independent observations, a known offset, an independent Gaussian nuisance covariate, and a single Gaussian kernel predictor. The design and generating response match the earlier 12-dataset pilot. The outer dataset and inner bootstrap seeds are new, fixed, and distinct. All datasets enter the reported denominators.

| Scenario | True scale (m) | Kernel coefficient | Generating distribution |
|---|---:|---:|---|
| Interior | 240 | 0.95 | NB2, theta = 2 |
| Weak effect | 240 | 0.12 | NB2, theta = 2 |
| Near upper bound | 730 | 0.95 | Poisson limit of NB2 |

The fit searches scales from 50 to 800 m. It uses the `joint_nb2` backend with the exact Poisson face. The optimizer uses the four declared pilot starts at 75, 200, 430, and 710 m and retains its existing numerical checks. The inference call uses `n_starts = 2` for its extra starts and retains the other starts used by the routine. No packages or model code are changed by the benchmark.

| Mode | Datasets per scenario | Bootstrap pairs per dataset | Total planned datasets | Total planned pairs |
|---|---:|---:|---:|---:|
| `bounded` | 8 | 49 | 24 | 1,176 |
| `precision` | 48 | 99 | 144 | 14,256 |

Only the bounded mode was requested for execution. The pilot took about 4.4 minutes for 228 pairs with one R process. Scaling by the number of pairs puts the bounded study near 23 minutes, with actual time depending on the optimization paths. The precision mode is a separate longer study with different dataset seeds. At pilot throughput it would take roughly 4.6 hours. Neither design changes the 90% target.

## Decisions and uncertainty

For a completed dataset, the test p-value is `(1 + b) / (B + 1)`, where `b` is the number of finite paired likelihood-ratio statistics at least as large as the observed statistic within the package's 1e-8 tie tolerance. At level 0.90, p-values greater than 0.10 include the scale. Equality rejects. The runner independently checks the p-value, exact integer decision, and adjusted 95% Monte Carlo limits against the package output.

Any failed pair makes the pointwise result unavailable. A case is Monte Carlo ambiguous when its adjusted lower limit is at or below 0.10 and its upper limit is above 0.10. Ambiguous cases retain their nominal decisions. The report also gives the sensitivity range from treating all ambiguous and unavailable cases as either rejecting or including, with clear decisions fixed. This range describes the observed cases.

Within each scenario, `nominal_rate` uses all datasets with a pointwise result. `returned_acceptance_rate` is the fraction of all attempted datasets that returned acceptance. Both include 95% exact Clopper–Pearson intervals, as implemented by [R's binomial test](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/binom.test.html). Those intervals reflect outer simulation uncertainty. The p-value limits reflect the inner bootstrap draws.

## Reproduce

Run this command from the radius development worktree with the path to the corrected package source and its compiled DLL. Use a stable source copy if development is ongoing. The reported run used a snapshot of the integrated root worktree.

```text
Rscript tools/benchmarks/nb2_bootstrap_replicated/run.R SOURCE OUTPUT bounded
```

Use a different output directory for `precision`. The script uses one R process and no parallel fitting workers. Repeating a command resumes completed datasets after verifying the source and design fingerprints. For a controlled pause, place an empty `STOP` file in the output directory. The current dataset finishes, then the runner writes a summary with the remaining rows marked pending. Remove that file before resuming.

- `plan.csv` contains all planned outer and bootstrap seeds.
- `datasets.csv` records every planned dataset and its outcome, including failures, Monte Carlo ambiguity, and warnings.
- `pairs.csv` contains every attempted bootstrap pair and its fitted statistics, convergence indicators, and messages.
- `cases/` holds atomic, compact per-dataset RDS checkpoints with outer optimizer records, profile audits, and pair outcomes.
- `scenario_summary.csv` gives denominators, exact binomial limits, and observed runtime.
- `configuration.rds`, `source_sha256.csv`, and `session.txt` preserve the source, design, and environment identifiers.
