# Probability precision availability: bounded software batch

Prepared while the declared spatial reference panel runs. No package source
edit, R test or sampler launch is part of this preparation. Wait until all four
spatial outcomes are retained and that study is safely frozen before numerical
validation or source edits. Do not overlap or restart the active controller.

## Observed issue and intended distinction

Both new null references retain 76 own-median binary events with estimate 0.5,
finite positive mean MCSE/ESS and undefined indicator Rhat. The current report
conditions all precision availability on finite indicator Rhat as well. The
scientific probability-MCSE criterion and auxiliary indicator-Rhat availability
need separate treatment. Neither null run can qualify: their raw diagnostics
retain five and three divergences. Every original study output must remain.

Inspect the actual installed posterior implementation before attributing the
undefined Rhat to balanced ties or folded diagnostics. Record its version and
reproduce the distinction with small deterministic original-chain matrices.
No new sampler, scientific dataset, prior, coordinate transform, fitting target
or numerical acceptance cutoff is declared. This is a software/reporting fix.

## Proposed change, contingent on the independent reproduction

Keep diagnose_bayes_probabilities' public arguments and numeric event estimates,
MCSE, mean ESS, Rhat and per-chain probabilities. Separate finite positive mean
precision from indicator-Rhat availability. An undefined Rhat must stay NA and
receive an explicit availability/reason field, rather than silently becoming a
finite substitute. Explain any intentional correction to precision_available
and its optional MCSE target decision in help and NEWS. That decision remains
computational precision, never a new fit-passed flag.

Preserve constant indicators, nonfinite scales, undersized original chains and
unavailable sampler reports as unavailable. Keep the original raw and derived
diagnostics, all priors, model densities, fitting/prediction interfaces, starts,
limiting/cold checks and uncertainty distinctions byte-identical. No fallback
after failed convergence, no dropped iterations, no changed Rhat/ESS/divergence
threshold, no new dependency, no staging/commit/version bump/merge/push.

## Focused assertions and validation sequence

1. Reproduce the old report against the checked installed package with a
   balanced binary event. Compute posterior::mcse_mean, ess_mean and rhat on
   the independent ordinary indicator matrix, and inspect installed Rhat source.
   Retain the old availability failure and all source/library identities.
2. Add a portable regression test: numeric diagnostics match posterior directly;
   mean precision follows its own finite positive MCSE/ESS, while Rhat availability
   follows the actual installed result. Do not require every future posterior
   version to return NA for this case. Test raw failure retention and no overall
   passed flag, target decisions, RNG/object immutability and chain order.
3. Retain constant, stuck-chain, nonfinite, single-chain, too-short, missing/duplicate
   parameter-layout and scalar/bracketed-name tests. Undefined auxiliary Rhat
   cannot certify convergence or rule out unseen tails/modes.
4. Compare predecessor/current numeric output and public formals, then run focused
   prediction/scale/derived-diagnostic/export compatibility tests. Inspect exact
   source changes against the 311-file baseline; preserve all unrelated edits.
5. Generate Rd with roxygen, preserve NAMESPACE exactly, update the computation
   vignette's explanation and execute it. Use the existing fast compact example;
   do not fit a model for documentation. Local browser preview remains policy-
   blocked; do not bypass that restriction or claim visual inspection.
6. Build a fresh dependency-complete source archive, check examples/vignettes,
   and run all unchanged installed tests separately with NOT_CRAN=true under
   the resource/interpreter cap. Old checks cover frozen old sources only.

A separately labeled research sidecar may later apply the tested mean-precision
distinction consistently to all four stored spatial records using their retained
MCSE/ESS/constant/Rhat fields. Keep original flags and reference decisions intact.
Apply all unchanged raw/derived and probability-MCSE gates, then the predeclared
paired threshold screen. Never manufacture a qualified null pair or rerun a fit.

All numerical work requires known live identities, at most two owned fit workers
and four actual R interpreters globally, one thread each, RAM >12 GiB/disk >5 GiB.
Full checks require an uncontended foreign-R-free window. No long launch after
04:00 Eastern October 29; first wake >=04:15 stops launches, checkpoint by04:35,
delete the one-night heartbeat by04:45, and no numerical work beyond04:45.
