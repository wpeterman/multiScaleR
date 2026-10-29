# Next software investigation: binary Rhat versus probability MCSE

The completed null-a-v2 run retains five divergences and is unqualified under
unchanged raw gates. Its derived-scale checks pass. All 228 fixed-threshold
precision reports pass, but 76 own-median events are marked unavailable by the
new adapter: each has estimate exactly 0.5, finite positive mean MCSE/ESS and
undefined indicator Rhat. The stored mcse/ess/rhat/reason fields preserve this
distinction. No iteration or chain needs to be rerun to retain this observation.

The scientific protocol requires finite nonconstant event MCSE <=0.005 and the
original raw/derived Rhat gates. Indicator Rhat was requested as a report field;
it was not given an additional pass cutoff. The current sr_event helper and the
new public probability diagnostic also make nonfinite indicator Rhat invalidate
all precision availability. Investigate that software condition independently.
Do not assume an undefined auxiliary diagnostic means MCSE is undefined, or that
finite MCSE establishes convergence or discovery of an unvisited mode.

After the remaining three declared fits finish, retain all four original result
files and freeze the study before any package source edits. In a new bounded
software batch, inspect the installed posterior implementation and use small
deterministic original-chain matrices to reproduce exactly balanced binary
events. Check mean ESS/MCSE and rank/folded Rhat independently, alongside stuck,
constant, single-chain and nonfinite cases. No new sampler configuration, seed,
simulation panel, target or prior is authorized by this investigation.

If a correction is supported, distinguish mean-precision availability from
indicator-Rhat availability in the new diagnostic's report/help/tests. Preserve
all original sampled/derived gates and raw diagnostics. A corrected probability
MCSE target remains computational precision, not a new fit-passed flag. Apply
any research-summary correction consistently as a separately labeled sidecar
for all four fits; never overwrite original decisions or rerun completed fits.
The null pair cannot qualify while its first run retains five divergences.

Until validated, keep all original reference decisions false. The compact
consumer may describe finite MCSE with unavailable Rhat without changing any
qualification. Changes to package source require a fresh source identity and
appropriate focused/documentation/full-package validation; the 5,016-assertion
check covers its original 311 files only.
