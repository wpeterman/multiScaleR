# Replicated NB2 pointwise results

The independent point panel finished all 96 prespecified datasets and
9,504 paired bootstrap refits. There were no failed pairs or unavailable
tests. The 99-draw bootstrap tails, Monte Carlo limits, and decisions all
matched independent reconstruction. Eleven pairs retained warnings.

| Scenario | Bootstrap inclusion | Exact 95% interval | Chi-squared inclusion | Exact 95% interval | Paired disagreement |
|:--|--:|:--|--:|:--|:--|
| Interior | 24/32 | 0.566–0.885 | 24/32 | 0.566–0.885 | 0/32 |
| Weak effect | 27/32 | 0.672–0.947 | 26/32 | 0.636–0.928 | 1/32, bootstrap only |
| Near upper bound, Poisson limit | 27/32 | 0.672–0.947 | 29/32 | 0.750–0.980 | 2/32, chi-squared only |

Both tests made the same decision for every interior dataset. The exact
interval for inclusion in that scenario lies below 0.90. The two tests
differed in three datasets across the other scenarios. The exact
intervals for disagreement were 0–0.109 in the interior, 0.001–0.162 in
the weak-effect scenario, and 0.008–0.208 near the upper bound.

Monte Carlo limits crossed the rejection threshold in 4 interior, 5
weak-effect, and 5 near-bound cases. In the interior scenario all four
ambiguous cases were nominal rejections. If every ambiguous decision
could flip, inclusion in the observed interior datasets ranges from
24/32 to 28/32. The corresponding ranges are 26/32 to 31/32 for the
weak scenario and 25/32 to 30/32 near the bound. These are descriptive
sensitivity ranges for the observed datasets.

Using only the first 49 draws changed one weak-effect decision from
inclusion to exclusion when the remaining 50 draws were added.
The 49-draw results were 24/32, 28/32, and 27/32 in the three scenarios.
The comparison uses common datasets and draws. It does not estimate
what a fresh 49-draw run would do on those particular datasets.

The panel ran from 00:45:54 to 02:28:57 ET. These findings
support retaining experimental calibration language for the pointwise
test. They concern acceptance of the known generating scale at one
candidate. Coverage of returned finite-grid sets is addressed by the
separately seeded grid panel.
