# Resume only the three pending declared fits

null-a-v2 completed once. It retains five divergences, zero depth hits, passing
derived-scale gates and all 228 fixed-threshold precision checks. Its reference
decision remains false. The original driver additionally marked 76 own-median
CDF events unavailable: their binary estimates are exactly 0.5 with finite MCSE
and ESS but undefined indicator Rhat. Preserve those original outputs and study
this diagnostic distinction separately. Do not rerun null-a or select another
seed, change the model, weaken raw/derived gates, or replace its result.

The v2 supervisor failed before launching null-b because File.Replace could not
remove the existing checkpoint destination. Its checkpoint records the failure;
no R worker remains. The destination has cloud-file attributes, but the specific
cause of the replacement error is not established. Preserve the entire v2
producer and checkpoint. No cleanup, attribute change or sync-process action.

The continuation uses unique, numbered atomic checkpoint files, requiring no
replacement of existing destinations. It writes each transition exactly once
and retains all earlier transitions. Before launch, test two distinct checkpoint
destinations and freeze the continuation source and original completed-result
identities in CONTINUATION_BASELINE.json.

Only null-b, smooth-a and smooth-b may launch, in their original order, using
the unchanged run-fit-v2.R, original seeds/settings/starts/priors/gates and -v2
output names. Require every pending output/process/exit path to be absent.
Reject any attempt to replay the completed first fit. Include the original
failed pre-sampling 21.2220588 seconds and completed null-a-v2 796.8090446 seconds
in the same cumulative 10,800-second wall budget. All remaining per-phase,
resource, exact-identity, one-worker/serial-chain and overnight limits persist.

Any later diagnostic software correction must be tested independently, apply
consistently to all four planned slots, retain original flags and separate
indicator-Rhat availability from mean-MCSE availability. It cannot erase five
divergences or qualify the null pair under the unchanged two-run requirement.
No further sampler configuration or scientific case is declared here.
