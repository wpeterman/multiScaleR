# Larger spatial response-field numerical study

This study compares the package calculation with a separately normalized
posterior simulation reference at known parameters. It does not estimate
model parameters from simulated responses. The joint target is the density
or mass of an entire assessment block conditional on the analysis outcomes.
The marginal targets are response means and total observation variances.

## Prespecified configurations

The main design includes 24 count scenarios. Poisson and NB2 cross two latent
SDs (0.5 and 1.8), three layouts and two independent outcome replicates. The
first layout uses eight analysis nodes and four assessment nodes with one
visit and nearest-neighbor correlation 0.25. The second and third use 16
analysis and six assessment nodes with nearest-neighbor correlation 0.75 and
one or four visits, respectively. Replicates use exponential and Gaussian
covariances. Four supplementary Bernoulli and Gamma cases use the latter two
layouts at latent SD 1.8. NB2 size is 3 and Gamma shape is 4. All seeds are
recorded. Four previously identified difficult quadrature cases serve as
independent positive controls for the simulation reference, rather than
as additional large fields.

## Independent reference and diagnostics

Sample the latent Gaussian posterior by elliptical slice updates at
assessment likelihood powers 0, 0.25, 0.5, 0.75 and 1. At each
temperature, independently optimize a numerical Gaussian reference
proposal. Elliptical slice sampling corrects for its density in the target.
If that proposal cannot be constructed, record its error and use the original
prior as the proposal for the same target. Each temperature has three
independent chains, 350 warmup updates and 1,600 retained updates per chain.
The temperature-zero posterior only fits the analysis responses. New field
innovations are sampled from their exact conditional Gaussian distribution
for the joint score. Marginal moments integrate those innovations
analytically or by univariate Gaussian quadrature for Bernoulli.

Bridge adjacent temperatures using the observed log likelihoods and an
equal-sample bridge. If the left and right posterior likelihood powers differ
by a, the normalizer ratio r solves the empirical equality between the left
mean of logistic(a times the log likelihood minus log r) and the right mean
of logistic(log r minus a times the log likelihood). The negative sum of
the adjacent log ratios is the joint block loss. Estimate Monte Carlo precision from 100-draw batches and
from independent chains. Use the larger standard error. The independent
chain and batch calculations preserve dependence between neighboring bridge
estimates. A score is resolved for this study if all split rank-normal R-hat
values are at most 1.05, all adjacent mean bridge overlaps reach 0.02, and its
standard error is at most 0.015 nats per assessment observation. A response
mean is resolved when the analysis posterior passes that R-hat criterion and
its standard error is at most 3% of the mean or 0.01 probability units for
Bernoulli. Predictive variance is resolved at a relative standard error of 10%.
All unresolved, failed and completed cases remain in the output. The
thresholds define how the reference is summarized in this study. They do not
change the package or define acceptable scientific errors.

The reference targets the exact model posterior at the frozen parameters,
with finite Monte Carlo precision. It does not verify a fitted model,
estimation uncertainty, data transfer, or calibration on ecological data.

## Prespecified precision follow-up

Rerun every larger-field case with an unresolved joint score or moment
reference, using independently seeded chains. For unresolved joint scores,
use nine likelihood powers at intervals of 0.125, with 3,200 retained
draws per chain. For unresolved moments, sample the analysis posterior
with 6,400 retained draws per chain. Preserve primary and follow-up
references separately. Also rerun the NB2 case with the largest resolved
joint discrepancy as a targeted replication of that tail result.
