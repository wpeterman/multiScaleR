# Shared prediction contracts

October 28, 2026, codex/inference-dev. The user authorized continued integration
through 04:45 Eastern October 29. This prediction-contract batch is validated;
the package check and complete installed tests ran in separate phases.
No staging, commit, version bump, merge or push.

## Implemented behavior

R/prediction_contract.R supplies one descriptive schema and small posterior,
analytic and native adapters. Bayesian prediction/scale summary tables and
posterior predictive simulations carry a prediction_contract attribute.
Prepared analytic/native snapshots retain a named prediction_contract list for
their prepared output types. Private factorized reconstruction records that its
quantiles condition on one supplied population vector.

The metadata distinguishes expected responses, linear predictors, kernel scales,
replicate responses, conditional choice probabilities and native components.
It records known/new/zeroed effect policies, joint-vector dependence, uncertainty
conditioning, omitted uncertainty, original-parameter block dimensions, retained
diagnostic status and lack of general calibration. New-group scale/mean quantiles
remain distinct from observation replicates. Zeroing is never labeled integration.
Native state/detection/initial-occupancy/transition meanings remain native.

Metadata cannot establish numerical eligibility. Analytic joint-covariance
availability is recorded separately and never bypasses the original checks.
The structural parameter-map flag does not validate coordinate transformations
or reconstruct integrated Gaussian coefficients. Call-level metadata can remain
on a subsetted data frame; the original row-level flags describe individual rows.

Existing public formals, priors, model densities, optimizer starts, diagnostic
thresholds, numerical predictions and random-number calculations are preserved.
Raw posterior draw outputs and cached native values/attributes are unchanged.
The intended visible addition is descriptive metadata on summaries, simulations
and prepared objects, documented in help and NEWS. The snapshot byte budget is
checked after metadata is attached. No new dependency or inference front end.

## Completed validation

* unit02:108focused contract assertions,0failures/errors/warnings/skips.
* Eleven existing prediction-related test files:571assertions,0failures/errors/
  warnings/skips using their latest successful per-file results. Nine files
  passed in compatibility01; the joint-mixed and conditional-logit modules
  passed their compatibility02 recheck. Both runs have identical full source
  maps. These are not two complete independent passing suites.
* parity01:8family/kernel cases across Gaussian/exponential kernels and Gaussian,
  Poisson, binomial and NB2 responses. Fitted/new-row means, link values, raw
  draws, scale summaries/draws and predictive replicates match exact retained
  predecessor values after removing only the new metadata. Seeded RNG/object
  immutability,6public formal lists and private fixed-vector reconstruction also
  match. No sampler study or calibration claim.
* documentation01:roxygen regenerated5help topics; all pass checkRd. NAMESPACE
  is byte-identical. The revised prediction vignette renders from current code
  using compact stored example fits; HTML is73,587bytes. No new fitting required.
* Review patch:13existing files changed and3new files, with exact predecessors
  retained for each existing file. All unrelated baseline files remain unchanged
  at the pre-check audit. The check freezes a complete package source rather
  than selecting a subset of the dirty tree.

The vignette-enabled package check and full installed tests completed as
check04/installed04, with the308-file
selected source manifest in check01/ and byte-identical archive reuse recorded
in check04-input.json. The standard check excludes its nested test subprocess;
all unchanged archive tests then run directly against its installed package.
The standard package check reports0errors,0warnings,0notes, with rebuilt vignettes
and examples. The initial installed suite passed4,678assertions with12existing
warnings and17skips. Seven skips were CRAN-only blocks, so supplement05 ran the
four affected unchanged archive files against the exact same installed package
with NOT_CRAN=true. All320supplemental assertions pass without warnings/skips;
193are previously uncovered. Final distinct coverage is **4,871passing assertions,
0failures,12existing test warnings,10development-only skips**. Do not add all320
or the focused suites to that total.

The read-only final audit passes all308source identities, additional-file scan,
exact archive identity, preserved prior evidence and empty index. All46sampled
resource inventories are complete with no foreign R hosts; observed maximum is
3actual interpreters, minimum17.12GiB available RAM and51.84GiB free disk. All
owned numerical work has exited. final-source-evidence.json freezes this scope
before the next batch; subsequent source changes require their own validation.
The old unified checkpoint correctly becomes stale after this deliberate batch;
it remains intact and must not be silently refreshed.

## Retained failed attempts and operational limits

unit01's105preliminary assertions passed before final refinements, but inherited
C.UTF-8 locale settings emitted startup warnings. Later workers clear those
settings only in their own process environment.

compatibility01 passed562assertions but had2conditional-logit test-runner errors.
The runner stored hashes in a global variable named identity, masking the base
function supplied as a test error handler. The package/test behavior was not
changed to hide those errors. The runner now uses source_hash_before and retains
the failed attempt. That first run also loaded a system glmmTMB built against
older TMB. compatibility02 uses the existing phase31 isolated glmmTMB1.1.14 with
TMB1.9.25 and RTMB2.0, matching the prior package-check setup, with no version
warning. No package installation was modified.

Browser visual inspection remains pending. The first preview connection lacked
its trusted bridge. The alternate browser action was rejected by automatic
approval review because local file URLs are not permitted; only HTTP/HTTPS are
allowed. No alternate route or policy workaround was attempted. Rendering,
generated help validation, and numerical example execution completed. The
rendered guide is retained at rendered/bayesian-prediction.html.

check01 was interrupted after its parallel tests exceeded the global interpreter
budget. The extra R supervisor around R CMD check caused avoidable interpreter
overhead. Its completed build/archive, examples, partial tests and shutdown
records remain retained; it is not relabeled as a successful full check. Only
exact root/descendant identities were stopped. Incomplete resource observations
are retained and do not establish uninterrupted isolation. check02 launches
R CMD check directly from PowerShell, reusing the exact built archive after a
fresh308-file identity check and R-free admission. The module audit counts loaded
R.dll interpreters separately from verified installed Rscript forwarders; the
four-interpreter limit is unchanged. check02 stopped fail-closed on an owned
R.exe CMD BATCH forwarder not yet classified by its monitor. Standard check
orchestration retains three R hosts before the two-worker tests, so check03
separates the check and test phases. The improved monitor recognizes a BATCH
forwarder only after verifying its installed binary hash, absence of R.dll and
its Rterm child. No parallel tests or gates were removed or weakened.
check03 exposed a separate environment NOTE for README/NEWS because native
R CMD check could not find Pandoc on PATH. Its interrupted record is retained.
check04 adds the already installed, verified Pandoc3.10 directory to its own
PATH; no global installation or package documentation rule is changed.

No performance or memory-gain claim is made for the new metadata. Prior snapshot
timings retain their original source and fixed-output scope. Study outputs are
compact JSON/logs, old source hunks and one small rendered guide; the full check
retains the dedicated source/check trees for each operational attempt. No chains, candidate populations,
samplers or tapes are serialized. No cleanup or historical blocked deletion.

## Overnight continuation

The bounded heartbeat integrate-multiscaler-until-4-45-am continues this chat,
not the three old chats. RESUME.md and NEXT_BATCH.md identify the next safe work.
This candidate's check/source audit is complete. The next small batch is
functional probability precision using retained unpermuted-chain diagnostics
and synthetic deterministic fixtures, without launching new NUTS studies.

Inspect live processes/resources before every numerical launch; never reuse a
historical PID as ownership proof. Max2owned workers/max4actual global R
interpreters including rsession/Rserve, one thread each, RAM>12GiB/disk>5GiB,
unknown inventory closes admission. Fullcheck admission was foreign-R-free;
check04/installed04-resources.jsonl record subsequent sampled ownership/resources.
No long launch after04:00; numerical launches close04:15; checkpoint by04:35 and
delete the heartbeat by04:45. Preserve all completed outputs and foreign jobs.
