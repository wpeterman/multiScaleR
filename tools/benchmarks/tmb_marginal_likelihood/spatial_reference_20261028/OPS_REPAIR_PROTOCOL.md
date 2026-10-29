# Pre-sampling driver repair

The first null-a execution failed before sampling. Its started/result/exit/log
records and the original frozen producer sources remain unchanged. The error
was "attempt to select less than one element in OneIndex". Source, input and
original producer identities were unchanged. No NUTS chain was produced.

Installed tmbstan 1.2.1 maps initial values using the names of RTMB's
last.par.best. The new driver's finite-value/gradient preflight evaluated
unnamed vectors, causing RTMB to drop those names. The repaired preflight keeps
the original block names, then passes exactly the same four numeric starts to
the sampler. A separate dry check reproduces the lost-name behavior, verifies
preserved names, all four identical numeric starts, finite values/gradients and
the installed tmbstan initialization indexing operation. It performs no sampling.

The supervisor also failed to replace an existing checkpoint because its original
Move-Item did not overwrite. The repaired writer uses adjacent temporary files
and File.Replace for existing destinations, tested before numerical admission.
Both errors and the original frozen scripts remain retained. No package source,
model density, prior, parameterization, sampler setting, seed, numerical gate,
functional threshold or planned scientific case changes.

The first standalone atomic-writer probe also failed before any R launch because
PowerShell converted a null backup-path argument to an empty string. Its source,
failure record and temporary replacement are retained. The corrected writer uses
NullString.Value; atomic-writer-validation02.json confirms repeated replacement.
repair-unit01.json confirms the named preflight and identical numeric starts.

The four planned scientific slots remain null-a, null-b, smooth-a and smooth-b.
The revised execution outputs use a -v2 suffix and explicitly retain the failed
pre-sampling null-a attempt as operational history. Exactly four NUTS fits remain
planned; none has run yet. Do not count a second driver execution as a new
dataset or discard the failed attempt. No concurrent or completed fit is replayed.
The cumulative wall budget includes the original failed phase's 21.2220588
seconds. Original 2,700-second per-phase and 10,800-second cumulative limits,
resource/identity admission and overnight deadlines remain in force.

Freeze the new repair producers separately in REPAIR_BASELINE.json before
launch, retaining BASELINE.json. Require both original and repair identity gates
before and after every phase. Preserve all failures, interruptions and pending
slots. No staging, commit, version bump, merge or push.
