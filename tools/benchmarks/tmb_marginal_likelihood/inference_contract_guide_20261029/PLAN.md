# Shared inference and prediction guide

Declared before drafting or running its new example. Start after the private
spatial source report's source/evidence batch is closed. This is one educational
documentation batch on codex/inference-dev, not a new fitting method or panel.

## Purpose and scope

Add one concise vignette, inference-contracts.Rmd, that connects the existing
optimization, Bayesian and private marginal-likelihood development tracks by
their targets and uncertainty contracts. A particularly important distinction
is that multiScale_bayes(method="marginal") still targets a posterior with its
declared proper population priors. It is not the private unpenalized marginal
likelihood or an AIC-ready likelihood. Laplace is an approximation technique,
not an inferential target. Existing conditional likelihood and native components
keep their own conditioning and interpretation.

Cover expected responses versus observation replicates; known/new/zeroed group
effects; joint Bayesian vectors versus conditional native standard errors;
analytic joint-covariance availability and its gates; snapshot reuse at fixed
rows and controls; functional MCSE versus convergence/importance diagnostics;
and kernel scale, scale-group SD and response-group SD units. Explain descriptive
prediction_contract metadata without implying eligibility or calibration. Public
spatial/Bayesian support is unchanged. Keep private spatial objectives and SD
traces clearly separate from supported public fits and confidence procedures.

Use existing source/help as the authority and link the quickstart, Bayesian
scales/computation/prediction and relevant public help topics. Avoid a new API,
dependencies, changes to R/src/tests or generated help, and a history of rejected
development experiments in user-facing prose. Add a brief NEWS entry and small
cross-links from existing relevant vignettes where useful. No invented citations,
em dashes or changes to numerical examples in existing guides.

## One bounded educational example

Use one deterministic small projected raster and six explicit point/response rows
to demonstrate kernel_prep() followed by multiScale_bayes(method="prepare") and
inspection of inference_metadata. No optimizer, chain, importance candidates,
scientific replicate, calibration dataset or speed claim is declared. Preparation
metadata must say posterior/proper priors and not_computed integration/uncertainty;
it cannot be interpreted as a fitted posterior. State natural distance versus
log-distance and observation-noise distinctions in the accompanying prose.

Run the example during vignette rendering with assertions about those prepared
metadata fields. Check help and vignette links against the current package,
restore any changed global state and keep the example small and self-contained.
No representative numerical posterior output is fabricated. Other worked fits
remain in their existing vignettes rather than being duplicated here.

## Validation and stopping rules

Freeze a baseline from the just-closed source and audit all predecessor evidence.
Render/execute the new vignette and inspect generated HTML content and runnable
assertions. Browser visual preview remains policy-blocked; make no visual-review
claim or alternate bypass. Generate documentation only if roxygen source changes
are actually needed (none are planned). Build a fresh candidate and run the
standard vignette/example check in an R-free window under the interpreter cap.
Verify R/src/tests/DESCRIPTION/NAMESPACE identities remain exact. The immediately
preceding complete installed-test evidence remains scoped to those identical
code/test sources; do not pretend it was rerun against the new archive. A fresh
full installed suite is required if any numerical/test/API source changes or a
new concern warrants it, otherwise avoid repeating numerical work for prose.

Known live inventory, one thread, max2owned fit workers/max4actual R globally,
RAM>12GiB/disk>5GiB; unknown inventory blocks launch. Focused phases max600seconds.
No long launch after04:00 Eastern, stop launches at first wake>=04:15, checkpoint
by04:35 and delete heartbeat by04:45. No staging, commit, version bump, merge,
push, broad cleanup, producer replay or additional sampler/simulation panel.

## Refinement before drafting or executing the new example, at02:03 Eastern

Source inspection found an existing bundled bayesian-annual-example.rds with
public preparation code and actual compact stored fits already used by the help
and vignettes. Reuse that example as a read-only consumer instead of constructing
the proposed six-point raster/response fixture above. The six-point proposal is
superseded and will not be run. This keeps one existing data/fit context across
the guides and avoids another example dataset or numerical fit.

Show method="prepare" with the bundled inputs/data and inspect the actual target,
proper-prior and not_computed uncertainty metadata. Then read the stored pooled
fit's expected-response prediction summary and its prediction_contract, using
only retained joint draws. Show actual metadata/output and explain omitted
observation noise and unchanged original diagnostics. Do not refit, sample new
chains, resample a new scientific cohort or generate new predictive replicates.
Simulation/functional precision concepts can link to the existing worked guides.
No new input data, fitted result, interval or numerical benchmark is invented.
