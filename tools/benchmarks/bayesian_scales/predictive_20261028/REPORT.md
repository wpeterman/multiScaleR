# Predictive calibration and numerical reliability

October 28, 2026. The planned fits and diagnostic investigations are complete. Nothing was staged, committed, pushed, merged, or versioned. Priors and diagnostic gates were preserved. This report distinguishes the frozen numerical study from concurrent changes to the working tree.

## Findings that affect use

An eligible fit passes all declared numerical diagnostics. A CDF gives the probability a model places below a specified value. Monte Carlo standard error (MCSE) measures numerical sampling uncertainty, which differs from calibration across independently generated datasets. AD denotes automatic differentiation, used to evaluate the posterior and its derivatives.

NUTS remains the better-supported method for prediction in these hierarchical models. It produced eligible fits on 54 of 56 primary inputs. The importance-corrected marginal approximation produced 10 eligible fits, 39 numerical rejections, and seven execution failures. All ten eligible marginal fits were Gaussian. A converged marginal grid alone did not establish a usable corrected posterior.

On those ten common eligible datasets, none of 190 usable predictive probability comparisons exceeded five combined Monte Carlo standard errors. Their new-year predictive widths and median recovery were close to NUTS. This is encouraging conditional agreement on a small, selected set. It does not establish equivalence, reliable calibration across response families, or safe use of rejected fits.

A numerical repair now allows all eight Poisson marginal fits to complete. Previously only one completed. All eight still fail the original importance diagnostics. The repair improves execution reliability without making unsupported accuracy claims. It is in the working source with regression tests and generated help.

## Design and provenance

The primary study reuses exactly 48 archived inputs and their original NUTS states from `prediction_20261027`. These include the original Gaussian/NB2 panel, the Poisson/binomial family panel, and an exploratory Gaussian follow-up. The cohorts remain separate throughout the analysis.

Eight independent prospective inputs add Gaussian and NB2 responses with an exponential kernel, four datasets per family. Each has 49 unbalanced training observations across five years, eight held-out observations, crossed site/year response intercepts, and an annual scale model with the fixed trend `time = (year_number - 3)/2`. The scale model is `~ 1 + time + (1 | year)`. Held-out tasks cross known/new year with known/new site. Inputs, independent truth calculations, ordinary densities, AD values, and fingerprints were audited before fitting.

This required 64 fresh primary fits: 56 marginal fits and eight NUTS fits. The 48 archived NUTS results were reused by reference. Two separately labeled NUTS diagnostic replications and 14 later stability/seed investigations bring the total to 80 compact new fit records. No diagnostic replicate replaces a primary result.

Both algorithms target the same proper posterior. The prior dictionary remains scale median 250, scale-intercept SD 0.6, scale-coefficient SD 0.5, scale-group SD 0.5, response-coefficient SD 1, response-group SD 1, and log-dispersion mean `log(2)` with SD 0.7. The exact dictionary is retained in the inputs and protocol. The fixed trend uses calendar labels rather than a newly estimated centering rule.

The marginal method uses importance-corrected marginal Laplace with nested response integration, `conditional_beta = FALSE`, t proposal df 5, minimum 65,536 and maximum 131,072 production candidates, and two refinements. NUTS uses four serial chains, 8,000 iterations per chain, 4,000 warmup, target acceptance 0.99, and depth 12. Each fit has the fixed 1,800-second internal budget. These are descriptive computational comparisons, not budgets matched for effective precision.

## Availability, including every failure

| Cohort | Response / kernel | Inputs | Eligible NUTS | Eligible marginal | Marginal rejected | Marginal execution failures |
|---|---|---:|---:|---:|---:|---:|
| Original | Gaussian / Gaussian | 8 | 7 | 2 | 6 | 0 |
| Original | NB2 / Gaussian | 8 | 8 | 0 | 8 | 0 |
| Family panel | Binomial / Gaussian | 4 | 4 | 0 | 4 | 0 |
| Family panel | Binomial / exponential | 4 | 4 | 0 | 4 | 0 |
| Family panel | Poisson / Gaussian | 4 | 4 | 0 | 0 | 4 |
| Family panel | Poisson / exponential | 4 | 3 | 0 | 1 | 3 |
| Gaussian follow-up | Gaussian / Gaussian | 16 | 16 | 6 | 10 | 0 |
| New trend panel | Gaussian / exponential | 4 | 4 | 2 | 2 | 0 |
| New trend panel | NB2 / exponential | 4 | 4 | 0 | 4 | 0 |
| Total | | 56 | 54 | 10 | 39 | 7 |

The two original unavailable NUTS references remain unavailable: `gaussian-004` was externally interrupted, and `poisson-exp-002` exceeded its time limit. All eight new prospective NUTS fits passed every fixed raw and derived diagnostic.

NUTS eligibility still requires finite Rhat at most 1.01, finite bulk/tail ESS at least 400, no divergences or depth hits, and all four finite E-BFMI values at least 0.3. Corrected marginal eligibility still requires Pareto k at most 0.7, raw ESS at least 1,000, indicator-tail ESS at least 400, maximum weight at most 0.01, and split-CDF discrepancy at most 0.05. Undefined diagnostics fail.

Among the 49 primary marginal fits that computed scores, 38 failed indicator-tail ESS, 34 failed split-CDF stability, 32 failed maximum weight, 21 failed raw ESS, and 12 failed Pareto k. These counts overlap. The dominant unresolved problem is weight concentration and tail precision, not merely grid convergence.

## Predictive recovery and calibration

Prediction uses complete joint parameter vectors. Each new year or site effect is sampled once per group and vector, then shared across its observations. Known effects remain attached to the correct rows. The coherent eight-observation total therefore retains dependence. Individual observation CDFs integrate response noise. New-year scale CDFs integrate the new annual effect conditional on its uncertain SD and use the entire fixed scale design, including the trend.

The tables retain 50%, 80%, and 95% intervals for expected responses, future observations, log scales, and the joint total. Population-SD intervals use their actual retained 90% endpoints. Dataset-cluster summaries, fixed-seed cluster bootstrap intervals, probability scores and computational precision are available in `summary.csv`, `paired-summary.csv`, and `report-tables.txt`.

The three final figures were rendered from completed results and visually inspected. They show [availability](figures/availability.png), [own-cohort coverage](figures/coverage.png), and [paired probabilities](figures/paired-cdf.png). `FIGURE_REVIEW.md` records the inspection and interpretation limits.

Selected 95% results for a new year and a new site are below. These are each method's own eligible cohorts, so the columns are not a paired method comparison.

| Cohort / response | Eligible NUTS / marginal datasets | Expected-response coverage, NUTS / marginal | Future-observation coverage, NUTS / marginal |
|---|---:|---:|---:|
| Original Gaussian | 7 / 2 | 100% / 100% | 100% / 100% |
| Gaussian follow-up | 16 / 6 | 100% / 100% | 93.75% / 100% |
| New Gaussian trend | 4 / 2 | 87.5% / 100% | 87.5% / 75% |
| New NB2 trend | 4 / 0 | 100% / unavailable | 100% / unavailable |

Visits within a dataset are dependent. For example, 87.5% across the four new Gaussian datasets has a conservative whole-dataset 95% bound of about 19.4% to 100%. Coverage of 75% in the two eligible marginal datasets has a bound of about 1.3% to 100%. Such wide bounds rule out a strong calibration claim. Passing 100% in a tiny panel is not evidence of ideal 95% calibration.

The ten common eligible datasets give the fairer algorithm comparison:

| Cohort | Common datasets | Marginal/NUTS 95% new-year log-scale width | New-year/new-site expected-response width | Future-observation width | Expected-response RMSE ratio |
|---|---:|---:|---:|---:|---:|
| Original Gaussian | 2 | 0.983 | 0.992 | 0.999 | 1.024 |
| Gaussian follow-up | 6 | 0.990 | 1.007 | 1.000 | 1.006 |
| New Gaussian trend | 2 | 0.981 | 0.985 | 0.998 | 0.993 |

These width ratios compare dataset-averaged widths. Their corresponding 95% coverage differences were zero. Recovery ratios this close to one, particularly with two datasets, cannot establish superiority. The 190 usable probability checks cover 19 targets per eligible input at the 95% reporting level, rather than counting repeated interval levels as independent evidence.

All four new NUTS Gaussian inputs covered the true new-year log scale, with median 95% log-scale width 4.222. All four new NB2 inputs also covered it, with width 3.974. These broad log-distance intervals reflect real uncertainty. Their exponentiated endpoints describe distances. A log-width is not a distance in map units. For the new Gaussian NUTS panel, the 90% scale-population-SD interval width was 99.4% of its prior width. The corresponding NB2 ratio was 91.7%. Response site/year SD widths contracted more, to about 68%/75% and 55%/73%, respectively. Five years provide limited information about between-year scale variation.

Some count-scale recovery errors are large because the preserved prior predictive distribution includes extreme counts. Their values are retained rather than filtered or stabilized through a new prior. A finite empirical RMSE or mean does not establish a well-behaved theoretical moment. Probability scores and log-scale summaries are useful complements.

NUTS probability and quantile precision uses ordered chains. Marginal CDF precision is conditional on the sampled importance proposal and cannot diagnose all missed modes or tails. Importance quantile MCSE remains unavailable. Constant probability sequences have unavailable precision rather than an invented zero. Archived discrete PIT randomization uniforms were not saved, so discrete paired PIT differences are omitted for those inputs. New inputs preserve shared uniforms. None of these omissions removes a failed fit or changes its eligibility. Binomial expected-response targets are success probabilities, while future observations are success counts under the declared trials.

The 37 primary computed-but-rejected pairs also had no five-MCSE flags among 559 usable checks, despite absolute CDF differences reaching 0.107. Large or limited precision estimates can make a discrepancy unflagged. This is another reason that an unflagged probability comparison cannot override a failed importance diagnostic.

## Failure investigation and repair

The exact original-seed `poisson-exp-001` reproduction captured one finite parameter vector with response linear predictors ranging approximately from -1,690 to 847. The proposal log density was finite, but the AD full log posterior was NaN. The independently evaluated ordinary full log posterior was negative infinity, correctly representing zero posterior mass after the Poisson mean overflowed. Moderate log scales ruled out extreme kernel scales as the cause of this counterexample. Only one small vector was retained, not a candidate population.

The new private helper rechecks an invalid AD target with the ordinary complete posterior at that same parameter vector. It accepts a finite ordinary value or negative infinity only. An unresolved target still triggers the original invalid-ratio error. The helper never drops a candidate, clips a parameter, alters a weight manually, changes a prior, or grants optimizer/curvature eligibility. It applies only to importance target evaluation. Adaptation and final samples retain their original random-number order.

The separate, versioned stability study reran all eight Poisson inputs and two controls at their original seeds. All eight Poisson fits computed. Across the rechecks, 338 invalid AD evaluations were recovered as zero-density evaluations, none needed a finite-value recovery, and none remained unresolved. The counts include pilots and discarded adaptation stages as well as production evaluations. Both Gaussian and binomial controls, which needed no recovery, produced exactly identical score rows and adaptation histories. The previously computed Poisson control was also identical.

All eight repaired Poisson fits still fail the fixed gates. Their ESS ranges from about 8.7 to 902. Even the best remains below 1,000, and their tail/stability checks also fail. The seven primary execution failures are retained in the primary analysis. The repaired results occupy a separate phase and do not retroactively improve its availability table.

Seven repaired Poisson fits have usable NUTS references on exactly the same inputs. Their 70 usable probability comparisons contain one five-MCSE flag. For `poisson-gaussian-003`, the row-8 new-year/new-site expected-response CDF at its true value 25.318 is 0.96428 under the marginal approximation and 0.93113 under NUTS. The difference 0.03316 is about 5.31 combined MCSEs (0.00624). Its central 95% upper endpoint is 36.27 versus 110.79, and its median is 0.706 versus 0.567. This fit has ESS 87, tail ESS 29.9, maximum weight 0.0616, and split-CDF error 0.286. The discrepancy is consistent with the unreliability already identified by those gates. Importance quantile MCSE is unavailable, so the endpoint difference is descriptive rather than an endpoint significance test. Across the repaired Poisson comparisons, absolute CDF differences reach 0.168. The eighth case has no usable NUTS reference after both timeouts. See `stability/nuts-comparison.csv` and `stability/nuts-paired.csv`.

Four independent-seed baseline diagnostic repetitions were declared separately. The Gaussian and binomial cases passed on their diagnostic seeds, the NB2 case remained rejected, and the Poisson case again stopped at invalid ratios. This shows seed sensitivity in proposal quality. Selecting successful repetitions would bias the study, so none replaced a primary result.

The two eligible independent-seed repetitions have 29 usable matched NUTS probability checks and no five-MCSE flags. The rejected NB2 repetition reaches an absolute CDF difference of 0.162. These are diagnostic replications, not extra independent calibration datasets.

The independent NUTS `gaussian-004` diagnostic passed in 77.61 seconds including scoring, supporting an external interruption explanation for its archival absence. The independent `poisson-exp-002` diagnostic again reached the fixed 1,800-second budget without returning usable draws. Its cause remains a runtime limitation on this case. No posterior comparison is possible. Two timeouts alone do not identify a specific geometric defect.

## Computational cost

Median fit plus scoring times, shown separately in seconds, were:

| Cohort / response / kernel | Marginal fit / score | NUTS fit / score |
|---|---:|---:|
| Original Gaussian / Gaussian | 40.62 / 27.86 | 104.14 / 2.65 |
| Gaussian follow-up / Gaussian | 20.79 / 12.80 | 85.05 / 4.41 |
| Original NB2 / Gaussian | 42.89 / 19.57 | 225.96 / 2.90 |
| New Gaussian trend / exponential | 36.36 / 23.16 | 99.66 / 4.50 |
| New NB2 trend / exponential | 48.65 / 21.09 | 321.75 / 2.36 |

These are computed-fit medians, not eligible-only paired speedups. Archival NUTS and current marginal timings occurred on a shared host, and the earlier Gaussian follow-up included extra auditing. Weighted prediction scored far more candidate vectors than the 16,000 NUTS vectors, explaining part of its cost. Rejected output is not a usable cheaper fit. Timeouts are retained separately, including both 30-minute Poisson references. No timing supports replacing NUTS across the count-response panels.

## When Bayesian pooling helps

The hierarchy is useful when years have uneven information, annual scales should vary around a meaningful common model, and prediction must include uncertainty about between-year variation. A weak year can borrow information from other years while its annual departure remains uncertain. New years require a population distribution for those departures. Full posterior prediction propagates uncertainty in its SD, coefficients, annual scales, response effects, and observation noise together.

Pooling can be harmful when exchangeability or the chosen trend is wrong. A sudden regime shift, residual temporal correlation, changing measurement, or weak habitat contrast may need a different model. New annual effects here are conditionally independent around the fixed trend, not an autoregressive forecast. Uncertain trend coefficients can broaden predictions farther beyond observed years even though the annual-deviation distribution is exchangeable.

This experiment compares two computations of the same Bayesian hierarchy. It does not add an unpooled or frequentist likelihood comparator and therefore cannot estimate a new pooling advantage over those alternatives. The broad intervals and nearly prior-width annual-scale SD intervals show why plug-in SD prediction can omit important uncertainty, but they do not by themselves prove superior accuracy. Use eligible NUTS when full uncertainty is central. Treat marginal output as provisional unless every diagnostic passes and predictive functionals have adequate precision.

## Integration and validation boundaries

The integration plan reuses the completed shared preparation, normalized density, inference metadata, and capability catalogue. It preserves likelihood versus posterior targets and native conditional uncertainty. Gaussian coefficient integration must reconstruct conditional coefficients once per joint vector before it can support prediction. Factorized annual likelihood requires uncoupled groups, and crossed response effects can violate that requirement. See `INTEGRATION_PLAN.md` for the ordered batches and ownership table.

At the coordination refresh, the Marginal Likelihood chat was developing a shared proper Gaussian block and a private spatial-scale prototype. Integrated likelihood was refining strict point-to-raster footprint matching. Their source and processes were left untouched. Those extensions require their own prediction covariance, target and calibration gates. The current annual study does not validate them.

Focused validation passed 705 assertions against the matched preserved source and 714 against the live source tested at 1:03 PM Eastern, including 41 new self-contained target-evaluation assertions. An additional 43 benchmark assertions cover the real overflowing-mean vector, finite targets across families/kernels, RNG preservation, unchanged public formals and unchanged acceptance code. Preparation passed 19 adapter assertions, nine later trend/discrete-PIT assertions, and ten resource-guard fixtures. An independent R audit checked score endpoints, coverage arithmetic, whole-dataset exact bounds, probability flags, archival PIT omissions, raw acceptance decisions, and all 80 retained result identities. Generated help passed `checkRd()`.

After the concurrent changes, a final focused compatibility run finished at 1:43 PM Eastern with all 714 Bayesian assertions passing, no failures/warnings/skips, and all 131 recorded source/test/native identities unchanged before and after the run. A separate R reconstruction matched all 720 repaired/diagnostic paired score rows, the 70 usable repaired-Poisson probability checks, their one flag, and all 14 per-result comparison summaries. These receipts update focused compatibility without certifying the full combined package or calibrating the new spatial prototype.

Earlier test bootstrap, plotting, and mismatched-source test attempts remain in the logs. The nine metadata failures in the mismatched test run resulted from using newer tests against the old numerical snapshot. The matched suite was rerun with its own tests and passed. No numerical tolerance or test was removed.

The frozen study passed 139 primary source/input identity checks and 152 stability-source checks. All 594 closed-study artifacts and 456 retained prior source files were unchanged. Concurrent owners subsequently changed `R/scale_model_capabilities.R` and `R/bayes_scale_prediction.R`. Their observed hashes are recorded separately. The earlier 714 assertions do not certify this later combined tree. Another owner is already checking a frozen integration candidate. This task did not launch a competing full check. A dependency-complete vignette-enabled check of the eventual selected combined source remains an integration gate, with exact fingerprints and a genuine resource window.

At the 1:38 PM coordination observation, that owner's check reported one test failure concerning the expected export `prepare_unmarked_data`, with 4,605 passing assertions. Its vignette stages passed. This is an owner-coordinated integration gate, not a clean combined-package receipt. The checked candidate was frozen before this task's numerical repair was promoted, so it cannot certify that later repair either. No failing test was edited here.

## Storage and closeout

The 80 new result records total 486,356 bytes, with a maximum of 9,064 bytes per record. They contain no chains, sampler objects, AD tapes, or candidate populations. Existing inputs and old sources were reused by reference. A small isolated repair source and two private reader wheels support reproducibility without changing package dependencies.

Eleven verified, unloaded temporary DLL folders from this task were deleted after owner-exit, retained-copy hash, sole-content and absolute-path checks. This reclaimed 83,728,503 bytes, about 79.9 MiB. One other registered folder was already absent. Source copies, all results, failures, logs and unrelated temporary files were preserved. The cleanup journal gives exact paths and retained hashes. The final cleanup audit confirms all 12 registered temporary copies are absent.

Across 1,474 owned-work resource observations, the conservative actual-R count never exceeded four, available RAM stayed above 13.7 GiB, and free disk stayed above 53.3 GiB. The new study is about 5.4 MB, with a separately retained 17.2 MB repair source, 1.7 MB reader wheels and a 0.06 MB documentation project. Free disk at the closeout observation was about 53.6 GiB. These are scoped observations on a shared host, not an attribution of every disk or memory change to this task.

All owned numerical workers exited. No foreign process was stopped. Resource inventories and the final artifact manifest accompany the report. Later foreign checks exceeded this study's conservative interpreter-count admission limit, so additional owned R work was deferred rather than bypassing it. The next approximation improvement should target proposal tail/weight stability on a newly declared independent evaluation panel. The present failures must remain visible.
