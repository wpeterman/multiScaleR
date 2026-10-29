# Shared model preparation and inference metadata

October 28, 2026. Implementation, exact predecessor parity, focused regressions
and help generation passed. A new full package check remains pending a genuinely
foreign-R-free window. No staging, commit, push, merge or version bump.
The closed simulation studies and their failures remain unchanged.

## Shared preparation

`R/scale_model_prepare.R` now holds the common row matching, formula designs,
kernel-bin preparation, fixed reference standardization and parameter indexing.
It takes an explicit reference scale and attaches no population priors. Existing
formula and prediction-design helper names remain for compatibility; those
helpers do not select an inference target or attach priors.

The existing Bayesian preparer validates the same prior settings, calls this
neutral layer and attaches the same priors and starting values afterward.
The private hierarchical-likelihood preparer calls the neutral layer directly.
It retains its reference-scale GLM initialization, identification checks and
exact zero-SD faces. Its preparation avoids population-prior construction entirely.

The default reference rules remain separate and unchanged. Bayesian preparation
uses its declared prior median. Likelihood preparation uses its existing explicit
reference or geometric reference default. Shared code does not silently force
different model specifications to use the same reference value.

Existing prepared-specification fields, order, formulas, row identities, scales,
contrasts and starting values must match their retained predecessors exactly.
This includes named scalar SD prior settings, cell-mean formulas, crossed groups
and supported multiple kernels. Formula environments stay frozen for compact
saved objects. No new large data copy, sampler serialization or AD tape is saved.

## Additive metadata

New fits contain `inference_metadata`, appended after existing fields. Public
arguments, classes, existing inference strings and diagnostic fields are retained.
Older saved fits can lack the metadata and continue through existing methods.

| Route | Target | Population uncertainty | Annual summaries or intervals |
|---|---|---|---|
| Bayesian prepare | Proper posterior specified, computation pending | Not computed | No intervals |
| NUTS | Proper joint posterior | Joint sampling | Marginal posterior credible intervals |
| Joint Bayesian Laplace | Proper posterior approximation | Joint Gaussian approximation | Approximate posterior credible intervals |
| Marginal Bayesian Laplace | Proper posterior approximation | Scale-SD quadrature and conditional Gaussian integration | Approximate posterior credible intervals |
| Importance correction | Proper posterior | Joint importance integration using the existing proposal machinery | Weighted/resampled posterior summaries with existing importance diagnostics |
| Private Laplace likelihood | Unpenalized hierarchical likelihood | Population parameters optimized | Conditional log-scale modes, back-transformed for distance estimates; calibrated confidence intervals unavailable |
| Private factorized likelihood | Unpenalized hierarchical likelihood | Population parameters optimized | Conditional log-scale modes and mean log scales; population uncertainty omitted |

The schema records its version, target, objective label, population-prior policy,
retained response normalizers, numerical integration by parameter block and
solver. Its uncertainty block states the interval type, population-uncertainty
treatment, annual-summary target and unestablished calibration status. Integration
descriptions apply when the corresponding block is present. The existing zero-SD
face record still identifies a selected model with a latent block absent.

Metadata describes a computation. It cannot grant numerical eligibility or
establish identification, posterior accuracy or interval calibration. Existing
diagnostics remain authoritative, including competing optima and limiting-model
witnesses. A conditional annual summary at fitted population parameters does not
propagate uncertainty in those parameters. Hierarchical likelihood still includes
the declared latent-effect distributions; its population parameters have no
Bayesian prior penalty.

## Validation scope and current status

PLAN.md declares predecessor/new preparation comparisons across two kernels,
four response families, grouped binomial trials, offsets, repeated landscapes,
row reordering, fixed/random scale formulas and joint multi-kernel designs.
Checks also compare ordinary/AD target values, gradients, Hessians, initial
parameters and zero-SD reductions. A syntax-tree check removes only metadata
assignments and arguments, then requires every remaining fit statement and
diagnostic criterion to match the predecessor.

Fourteen paired Bayesian/likelihood preparation cases match exactly, including
prepared fields, their ordering, parameter starts and zero-SD reductions. Four
additional Bayesian cases cover multiple kernels, crossed grouping terms and
named scalar prior settings. All ordinary/AD target, gradient and Hessian
comparisons have zero recorded difference. Removing only metadata assignments
and constructor arguments leaves six fitted-code modules with identical parsed
statements, including their original diagnostic criteria.

All 101 new preparation/metadata assertions pass. The existing Bayesian,
Laplace, marginal, factorized and shared-response tests pass 1,325 assertions
across 23 files. The combined total is 1,426, with zero failures, errors, warnings
or skips. The final regression run finished at 07:12 Eastern. Roxygen generated
the updated multiScale_bayes.Rd help from source, and its Rd check passed.
No generated help was edited manually.

Initial resource admissions were held by three foreign RStudio sessions and
another package's validation. After the user closed idle sessions, the numerical
work ran under the original worker, RAM, disk and single-thread gates. No foreign
process was stopped by this task. Full package checking still requires a
genuinely foreign-R-free window; the remaining RStudio session and concurrent
validation keep that gate pending. Earlier check results are not a full check
of this new batch. No archived NUTS or six-case simulation fit was repeated.

Two validation failures are retained. The first worker lacked the existing
private RTMB research library in its library path and stopped before derivative
tapes were constructed. Its successor used the same existing library as the
closed studies, without installation. One new backward-compatibility assertion
then assumed observations were prediction-matrix rows; the existing contract
uses draws as rows. The corrected test compares predictions with and without
metadata and verifies the documented dimensions. Production code was unchanged
after the successful parity result. Preparation parity includes the existing
GLM initialization, but runs no hierarchical scale optimizer or sampler.

validation-summary.json verifies the current critical source hashes against the
successful parity and regression records. Only the twelve intended existing
source/test/help files differ from the task baseline. No unrelated baseline
file changed. The new private helpers and their tests are additive. Compact
validation artifacts total about 0.21 MiB; no sampler, draw array or AD tape was
saved. Previously blocked temporary deletions were not retried.

## Boundaries and next steps

This batch shares preparation for the Bayesian and private hierarchical-likelihood
routes. It does not replace specialized profile, conditional-logistic or other
existing likelihood preparation. Their separate contracts remain unchanged.
It adds no model family, random slope, temporal covariance or automatic method
selection. The factorized independence restrictions remain in force.

After validation, the next integration step is a dependency-complete development
tree with capability checks and these explicit target contracts. Public transport
reconstruction, population-SD confidence inference and broader independent
calibration remain separate scientific priorities. The exact Gaussian beta
research prototype remains experimental and is not exposed by this change.
