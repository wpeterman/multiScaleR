# Shared capability contract for scale inference

October 28, 2026, completed at 11:47 AM Eastern. No staging, commit, push, merge
or version bump. This batch
continues the completed shared response-density, model-preparation and inference
metadata work. It leaves the archived calibration studies and failures intact.

## What changed

R/scale_model_capabilities.R defines a private seven-route table and a shared
structural preflight. The same route catalogue now validates inference metadata,
so target/method combinations cannot drift between those two layers. Metadata
values, field ordering and the public fitting signatures remain unchanged.

Existing marginal-likelihood preparation, factorization checks and the marginal
posterior's one-scale-SD check use this layer. Public unsupported marginal models
stop before solver dispatch. Direct private joint latent Laplace calls now reject
unsupported multiple-kernel or crossed scale/response-term specifications before
optimization, rather than relying on callers to have used the private preparer.
Existing restriction messages are retained. For unsupported requests, the
structural error can now precede an optional-software or solver-control error.

Model support is deliberately separate from package availability, resource
admission, identifiability, numerical eligibility and uncertainty calibration.
The capability table never selects a fallback, changes a prior, grants an
accepted fit or supplies missing confidence intervals. Private likelihood fits
remain private. Shared checks return the prepared specification unchanged.

The public help explains why marginal posterior integration needs one total
scale SD rather than one kernel, and why a second random scale term creates a
second SD even when both terms use year. INTEGRATION_STATUS.md gives the model
support table and remaining integration priorities.

## Validation scope

parity-v2.json passes 16 exact predecessor/current prepared-model comparisons:
two kernels, four response families and presence/absence of a response grouping
term. Factorization decisions and messages match their predecessors in every
case. All seven inference metadata routes are exactly identical. An explicit
parsed-code gate removes only the declared preflight additions/substitutions;
all remaining numerical, prior, start and diagnostic statements are unchanged.
Parity runs existing reference-scale GLM initialization, but no hierarchical
scale optimizer or sampler.

The new capability tests cover multiple kernels, zero/one/two scale SDs, crossed
scale terms, Gaussian response groups nested within annual groups, response
effects crossing annual groups, aliases with different grouping-variable names,
shared scales, exact zero-SD faces and orphan response latent parameters. Mocked
solver entries prove unsupported requests stop before fitting. Existing focused
tests supply response/objective derivative checks and method-specific diagnostic
regressions. No archived NUTS or six-case scientific fit is repeated.

The final evidence passes 1,598 assertions across 26 files, with zero failures,
errors, warnings or skips. This includes 172 new capability assertions, 101
preparation/metadata assertions and 1,325 existing regressions. Eighteen completed
files were reused only after exact source/test identity checks and a parsed
proof that the sole helper correction removed two generic checks irrelevant to
their valid prepared models. The revised capability/preparation/metadata tests,
the failed refinement tests and the four missing files ran on the final source.
reuse-evidence.json records that restricted reuse; no completed numerical fit
was rerun merely to fill the queue. No full-tree claim is inferred from it.

Roxygen regenerated help from source, tools::checkRd() passed, and the final
validation verified all current R source hashes against its start. final-audit.json
confirms those hashes still match. Seven intended existing source/help files
changed, with one new helper and one new test file. Concurrent changes in native
prediction, prediction help/tests and the main guide were retained. Their own
validation remains separate. No generated Rd file was edited manually.

A full package check still requires a genuinely foreign-R-free window. The
final inventory contained one RStudio host and two foreign R interpreters.
They were untouched. About 15.6 GiB RAM and 56.6 GiB free disk remained; no owned
root or numerical worker remained. Study files were about 0.23 MiB at audit.
There is no new automation or queued run. A check of an earlier source tree is
not a full check of this integration tree.

RETAINED_FAILURES.md records the original collinear test fixture, the resource-only
pause caused by incomplete details for a short-lived foreign R interpreter, and
the three optimizer refinement errors caused by the new helper requiring family
and kernel fields in structural-only fixtures. The original helper/test/logs and
progress remain. Family/kernel validation stays in neutral preparation; the
original refinement test was retained unchanged and now passes all 39 assertions.
The final successful controller exited zero. Resource retries use fresh inventories
and retain the same worker, RAM, disk and single-thread limits.

## Scientific and development limits

This batch adds no random slope, temporal covariance, response family or
general scalar factorization. Full Bayesian joint routes still have broader
structural support than the private likelihood prototypes. Both hierarchical
likelihood and Bayes can pool annual effects. Bayes additionally integrates
population-SD uncertainty under explicit priors; a structural capability check
does not show that those priors help recovery or that intervals are calibrated.

The closed studies' accepted Gaussian/sbc-111/sbc-057 decisions and rejected
prior_shift-018/sbc-123/sbc-187 decisions remain unchanged. Public reconstruction
for transported posterior draws, boundary-aware likelihood SD uncertainty and
broader independent calibration remain separate substantive priorities.

No temporary deletion or prior policy-blocked cleanup is attempted. Only small
source predecessors, logs and atomic summaries are retained, with no serialized
AD tape, sampler, draw array, package-tree copy or duplicate simulation dataset.
