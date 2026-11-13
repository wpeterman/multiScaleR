# Bayesian approximation validation

The primary question is whether importance-corrected marginal Laplace produces
a usable approximation to the same posterior as NUTS, across the specified
families and designs. Two importance proposals are compared at the same
declared candidate budgets. The comparison retains every fit, including
numerical rejections and execution failures.

## Design

The simulator draws the full parameter vector independently from the model's
declared priors and generates training and future responses with base R. It
recomputes the landscape kernel and the training reference transformation
without calling the package's likelihood or simulation functions. The design
uses Gaussian, Poisson, grouped binomial, and NB2 responses. Each family has
two panels:

- Gaussian kernel with three annual scale effects and no response random effect.
- Exponential kernel with a fixed trend, four annual scale effects, repeated
  site response intercepts, and an unbalanced observation schedule.

The panels include future observations at a new year and a new site. Future
observations share each simulated group effect, so a joint predictive density
uses one common draw per group. Comparisons include training log-scale
probabilities at predeclared 150 m and 350 m thresholds, log scale and response
coefficient prior-predictive ranks, and the joint future response log score.

Pilot and qualification seeds are disjoint. The pilot runs one input per
family and panel to check the workflow and characterize failures. The
qualification protocol uses 48 independent inputs per family and panel.
Do not combine panel sizes or use the pilot to assess general calibration.
Priors, posterior evaluation, and the package's diagnostic cutoffs are held
fixed across methods. The reference uses four independent NUTS chains and must
pass the package's raw and derived checks.

## Assessment

Evaluate each family, panel, and proposal separately. Report the number of
eligible, rejected, failed, and unrun fits, with exact binomial intervals for
eligibility. A cell needs a one-sided 95% lower eligibility bound of at least
0.90 before considering wider use. Assess prior-predictive ranks and empirical
central interval coverage using independent dataset-level replicates. Use an
exact binomial interval around 90% coverage and flag systematic rank deviations.
Compare fitted-scale probabilities against NUTS only where both fits and both
Monte Carlo precision estimates are available. A discrepancy exceeding both
0.05 and five combined MCSEs requires investigation. Also inspect prediction
log-score differences and retain their whole-dataset distribution. The
comparisons do not overwrite fit eligibility or establish field-data
calibration. A claim is restricted to the family, design, priors and settings
represented by evidence. No label changes on the basis of the pilot.

## Run

From the package root, load a matching compiled package or use `pkgload` to
build it. The runner needs the existing RTMB, tmbstan, rstan, posterior and
digest dependencies. Supply an external output directory. An optional case or panel prefix limits execution to that cohort.
Suffixes `":odd"` and `":even"` select independent replicate shards. Start one
shard before the other creates records so both share the saved configuration.
Run with `assess` to summarize available records without fitting. Rerun the
same command to resume. A different source or protocol fingerprint
requires a new output directory. Output consists of compact atomic RDS records,
a fit summary, and a probability comparison table. Raw candidate matrices and
Stan objects are not saved.

```text
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/run.R pilot OUTPUT
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/run.R qualification OUTPUT gaussian-annual:odd
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/run.R qualification OUTPUT gaussian-annual:even
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/run.R qualification OUTPUT assess
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/assess.R OUTPUT
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/verify.R OUTPUT
```

Use `qualification OUTPUT all` to run all eight panels serially, or run
the other family and design prefixes using the same odd and even shard
syntax. The finished [fit archive](results/replicated/fit-records.zip)
contains the compact RDS records and configuration for independent
assessment and fingerprint checks.

[RESULTS.md](RESULTS.md) records the eight replicated panels and the
separate eight-dataset pilot. The two largest predictive-score differences
can be reassessed with quadrature over the new effects and independently
seeded NUTS fits.

```text
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/diagnose_tail.R qualification OUTPUT assess
Rscript tools/benchmarks/bayesian_scales/qualification_20261112/diagnose_poisson_tail.R qualification OUTPUT assess
```
