# Spatial response-field CV quadrature study

This study compares the non-Gaussian field CV calculation with independent
numerical integration at known, frozen parameters. The target is the joint
assessment negative log mass or density conditional on the analysis responses.
The study also compares assessment response means obtained using analysis
responses only.

## Design recorded before execution

The factorial crosses four response families, field SDs of 0.25, 1 and 2 on
the linear predictor, nearest-neighbor latent correlations of 0.20 and 0.85,
one or eight visits per node, and three layouts. The layouts have one analysis
and one assessment node, two analysis and one assessment node, or one analysis
and two assessment nodes. The second layout has analysis nodes on either side
of the held-out node. Eight outcome replicates in each cell use predetermined
seeds. Four use exponential correlation and four use Gaussian correlation.
The design comprises 1,152 non-Gaussian scenarios. It also includes 192
zero-field cases and 288 exact Gaussian checks. Fixed predictors vary by
observation. NB2 size is 3, Gamma shape is 4, and Gaussian noise SD is 0.65.

The unit of replication is one simulated analysis and assessment split. All
attempted scenarios and their errors are retained. Accuracy summaries use
cases with converged independent references and report reference failures and
approximation failures separately. Report the median, 90th, 95th, 99th and
maximum absolute joint loss errors, per-observation loss errors, and
response-mean errors. Report absolute probability error for Bernoulli and
relative error for log-link response means. Treat 0.01, 0.05 and 0.10 as
descriptive error bands in nats per observation, probability units, or
relative mean units as applicable. These bands do not impose a universal
tolerance for scientific or model-selection decisions.

## Reference

Calculate the analysis and combined log marginal likelihoods using normalized
base-R conditional densities and adaptive tensor Gauss-Hermite quadrature in
one to three dimensions. An independent numerical BFGS optimizer supplies the
proposal center and Hessian. The quadrature includes the full normalized
prior-to-proposal density ratio. It does not use the package's latent mode,
derivatives, Hessian, or Laplace normalization. Both the integration order and
the proposal spread are varied before classifying the reference as converged.
The refinement and spread checks must agree within 0.000005 for both log
marginal likelihoods and within 0.000005 times the larger of one and the
reference mean magnitude for each prediction. Quadrature uses up to 80 nodes
per dimension.
For response means, integrate over the analysis posterior quadrature and the
conditional Gaussian innovations at the assessment coordinates. Log-link
moments use the Gaussian moment formula, and Bernoulli probabilities use a
separate one-dimensional quadrature. Exact zero-field references use direct
observation densities. Exact Gaussian references use a separate
observation-level Gaussian calculation.

These comparisons address the numerical approximation at frozen parameters.
They do not test ecological scale or range recovery, parameter uncertainty,
field variance-boundary fitting, interval calibration, or transfer to
high-dimensional fields. The package's native fitting and leakage fixtures
are separate tests.

## Follow-up after the initial assessment

For the three log-link families, run a separate experiment on the same
scenarios and compare the first two conditional mean moments and total
observation variance. The independent reference integrates the exact
conditional mean moments over the analysis quadrature weights. Both log
moments must agree within 0.000005 across order refinement and proposal spread,
using up to 112 nodes per dimension. Run 30 alternating-order timing pairs
at the training estimates from the executable Poisson vignette. Time the
baseline score and prediction and the incremental warm-started tilted solves
separately. Keep all moment outcomes, quadrature diagnostics and timing pairs.
