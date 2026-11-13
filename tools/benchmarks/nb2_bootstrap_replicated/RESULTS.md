# NB2 pointwise bootstrap results

The generating scale was included in 21 of 24 datasets at the nominal 90%
level. All 1,176 paired bootstrap refits completed. The package decisions,
p-values, and Monte Carlo limits matched independent calculations from the
retained likelihood-ratio statistics in every dataset.

| Scenario | Included / attempted | 95% exact binomial interval | Monte Carlo ambiguous | Observed inclusion if ambiguous decisions switch |
|---|---:|---:|---:|---:|
| Interior | 8 / 8 | 0.631 to 1.000 | 1 | 0.875 to 1.000 |
| Weak effect | 6 / 8 | 0.349 to 0.968 | 2 | 0.750 to 1.000 |
| Near upper bound, Poisson limit | 7 / 8 | 0.473 to 0.997 | 2 | 0.750 to 1.000 |

The exact intervals describe the outer sampling uncertainty for inclusion
under this 49-draw test. All three nominal rejections had Monte Carlo
intervals crossing the rejection threshold. A weak-effect case at p = 0.10
correctly rejected at equality. Two nominal acceptances also had Monte
Carlo intervals crossing the threshold. The five decisions remain in the
nominal counts above.

| Case | Bootstrap p-value | 95% Monte Carlo limits | Nominal decision |
|---|---:|---:|---|
| Interior 07 | 0.16 | 0.078 to 0.287 | Include |
| Weak 01 | 0.10 | 0.042 to 0.212 | Reject |
| Weak 03 | 0.08 | 0.033 to 0.185 | Reject |
| Near bound 04 | 0.06 | 0.025 to 0.157 | Reject |
| Near bound 05 | 0.12 | 0.053 to 0.238 | Include |

Availability and localization differ. Seven of eight near-bound datasets
selected the exact Poisson face. Six had a bound-limited likelihood profile.
Five of eight weak-effect datasets had an unresolved profile, and two
had a bound-limited profile. An included generating scale can therefore
accompany limited information about its location.

Three bootstrap pairs retained warnings. Of 28,227 combined constrained
and unrestricted optimizer starts recorded across those pairs, 28,171
were eligible. Every pair still had eligible fits. Thirteen datasets
retained an initial, fitting, or inference warning. See the pair CSV for
the warning text and the RDS checkpoints for the point-fit and profile
audits.

## Runtime and evidence

The bounded run took 27m29s from its first dataset to the last completion
on one R process. Median time per dataset was 92.5 seconds for interior,
64.7 seconds for weak effect, and 42.1 seconds for near bound. The pilot's
4.4 minutes for 228 pairs suggested about 23 minutes for this run under a
linear extrapolation. Using the bounded run's observed time per pair
projects about 5.6 hours for the separate 144-dataset, 99-pair mode. That
mode was not executed.

The completed archive contains 24 independently seeded dataset
checkpoints and all 1,176 paired outcomes. The run was executed from the
corrected source snapshot, and all recorded source hashes were
reverified. [Scenario summary](bounded/scenario_summary.csv) ·
[Dataset outcomes](bounded/datasets.csv) · [Pair outcomes](bounded/pairs.csv) ·
[Reproduction instructions](README.md)

These results verify execution and decision arithmetic under the three
specified simulations. Eight datasets per scenario leave wide inclusion
intervals. The study supports retaining experimental language for
calibration. It evaluated a test at the generating scale, so it supplies
no empirical coverage estimate for a returned finite-grid confidence set.
