# Analytic scale assessment results

## Independent GLM likelihoods

All 96 `profile_glm` fits and all 24 `joint_nb2` fits completed. The
independent same-scale reference was available for all 120 fits. At the
fitted scale, Gaussian, Poisson, and grouped binomial likelihoods matched
the native GLM calculations within 3e-14 negative log-likelihood units.
`joint_nb2` matched an independent normalized NB2 likelihood within 3.0e-9
units, including the exact Poisson face. It selected that face for five of
the eight Poisson-limit datasets generated near the upper scale bound.

Three near-bound NB2 `profile_glm` fits had gaps above 0.0001 compared with
the reference that allowed the Poisson face. At the package's finite upper
dispersion limit of theta 100,000, independent native coefficient fits
matched within 1.5e-8 negative log-likelihood units. The largest discrepancy
of 0.000500 was the difference between that bounded NB2 likelihood and the
exact Poisson likelihood. The independent coarse grid did not find a
better objective within the same likelihood target. The [face check](check_poisson_face.R)
records these comparisons separately.

### Scale uncertainty

In the interior scenario, the nominal 90% log-scale Wald interval
contained the truth in 25 of 32 `profile_glm` fits across four families and
5 of 8 `joint_nb2` fits. These small cells support the warning about
undercoverage and provide limited family-specific precision.

In weak-effect fits, `profile_glm` returned 16 of 32 Wald intervals, of which
13 covered the truth. `joint_nb2` returned two of eight, both covering the
truth. Near the upper bound, `profile_glm` returned 19 of 32 intervals, of
which 18 covered the truth. `joint_nb2` returned no Wald intervals in the
eight near-bound fits. The conditional coverage among available intervals
does not describe all generated datasets.

## Joint mixed point estimation

The `joint_mixed` backend returned 53 of 54 fits. One NB2 fit generated
under the Poisson limit failed its joint Hessian and cold native audits.
For one other such dataset the joint fit returned, but the independent
native refit at its scale had a nonpositive Hessian. The remaining 52
same-scale checks agreed within 4.8e-7 negative log-likelihood units. None
of the independently evaluated grid points yielded a higher likelihood
than the accepted joint fit.

These results support the point estimator for the evaluated Poisson,
NB2, and zero-inflated Poisson random-intercept models. The experiment
does not establish scale interval coverage. Native prediction and
coefficient standard errors condition on the selected scale.

## NB2 pointwise bootstrap

The 12-dataset pilot completed 228 paired bootstrap refits. It exposed
a rejection-cutoff error. At confidence level 0.90, subtracting the
level from one rounded below the represented p-value 0.10 and incorrectly
accepted equality. Comparing on the confidence-level axis restores the
intended decision. The corrected rerun reproduced all 12 p-values and
changed only the decision at p = 0.10. Eleven of twelve generating scales
were accepted and two tests had Monte Carlo intervals crossing the cutoff.
The before-fix summary is retained for provenance.

The profile and bootstrap panel uses 19 bootstrap samples per dataset.
Its p-values have resolution 0.05 and its recovery rates have substantial
sampling uncertainty. It assesses execution, reporting, and the cutoff
decision. Calibration of bootstrap confidence sets remains experimental.

Each run retained compact RDS outcomes, including failed starts and
independent reference warnings, in its external output directory. The
summary tables in this directory contain all attempted datasets. See
[README.md](README.md) for the fixed scenarios and reproduction commands.
