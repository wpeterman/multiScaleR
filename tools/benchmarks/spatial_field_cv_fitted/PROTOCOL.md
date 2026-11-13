# Fitted spatial response-field predictive assessment

## Question

How often does the public field-CV workflow return an eligible fitted score
when a shared Gaussian spatial field generates the data, and how do those
plug-in predictions compare with predictions at known parameters?
Each simulation independently generates one field, its training outcomes
and one spatially structured assessment block. Analysis and assessment
never share a location. The unit of replication is the simulated field.

## Frozen design

The main panel has three response families, three designs and 40 simulations
per cell, for 360 fits. Each seed identifies a simulation independently
of execution order. The same data are used for the fitted and known-parameter
predictions within a case. Families use a Gaussian identity, Poisson log,
or NB2 log link. The Gaussian observation SD is 0.65, NB2 size is 3.5, and
the field and fixed coefficients vary by design only as stated below.

| Design | Field covariance | Range (map units) | Field SD | Analysis and assessment locations | Visits per location | Geometry |
|---|---|---:|---:|---|---:|---|
| Interpolation | Exponential | 1.6 | 0.65 | 28 and 8 | 3 | Interspersed held-out locations on a 6 by 6 grid |
| Sparse strong field | Exponential | 3.0 | 1.25 | 17 and 7 | 2 | Interspersed held-out locations on a 6 by 4 grid |
| Buffered transfer | Gaussian | 2.5 | 0.85 | 24 and 12 | 3 | Training and assessment bands separated by 2 excluded columns |

Coordinates are in projected Cartesian units. Log-link expected counts
use intercept 0.55, habitat coefficient 0.50 and a log exposure offset.
Gaussian means use intercept 0.80, habitat coefficient 0.55 and the same
offset. Covariates and exposures are sampled independently of the field.
These are model-correct comparisons for specified designs. They leave
covariance misspecification, ecological scale estimation and nonstationary
fields outside the panel.

## Targets and evaluation

Fit the original public `cross_validate_field()` at each analysis split,
with the generating covariance family and the documented fit controls and
eligibility criteria. Keep both native starts and all failures. Each run
uses glmmTMB 1.1.14 and TMB 1.9.21. No fit uses assessment responses.
Successful fitted values include all public predictions and the joint loss.
Report fit availability over all 40 attempted simulations in each design
and family, followed by conditional paired comparisons for the scored
simulations. Keep the true-parameter predictions available for cases whose
fitted prediction fails.

Use the same normalized Laplace scoring and tilted moment integration at
generating parameters for counts. For Gaussian responses, calculate the
known-parameter predictive distribution by independent dense
observation-level Gaussian conditioning. Compare paired fitted minus
known-parameter joint log loss per assessment observation, per-field
mean error, squared error, predictive variance, and the mean squared
standardized prediction residual. Calibration of the last statistic
corresponds to a value of one. That is a moment check rather than a
distributional check. Report RMSE relative to predictive SD through the
ratio of summed squared error to summed predictive variance within
simulation. Report paired confidence intervals from resampling simulated
fields, retaining location and repeated-visit dependence. Field-bootstrap
intervals are descriptive with 40 replicates per cell.

For Gaussian responses, use normal CDF probability integral transforms
and 90% predictive interval coverage at each assessment row. Average row
statistics within each field before pooling across fields. For counts, use
a predetermined subset of first and second replicates in each design
and family for independent elliptical slice references at known and
eligible fitted parameters. The sampler uses three chains and exact
normalized observation densities. Use its marginal CDF and moment
references to check the prediction calculation at those parameters.
PIT values from such a small subset are illustrative. Compare 32- and
64-node quadrature for the conditional innovations. A PIT reference is
resolved if its split rank-normal R-hat is at most 1.05, its MCSE at most
0.015, and the absolute refinement gap at most 0.002. For the joint score
also require every adjacent bridge overlap to reach 0.02 and MCSE at most
0.015 nats per row. For mean and variance require relative MCSE at most 3%
and 10%, respectively. Keep all attempted references.

Coefficient and covariance recovery are descriptive for eligible fits.
Record estimates of intercept, habitat coefficient, field SD, range and
Gaussian observation SD or NB2 size. Report medians, 10th and 90th
percentiles of fitted-to-generating variance-parameter ratios and the
fraction outside one third to three times the generating value.
This range describes practical instability. It does not establish that
an individual fit is invalid. Retain failed starts and eligible fits
whose two starts disagree.

Run each planned case once with its original gates. Retain every failure
and denominator. Independent sampler checks are a separate cohort.
