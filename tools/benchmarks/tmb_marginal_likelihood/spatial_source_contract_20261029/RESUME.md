# Spatial source preflight: closed and validated

Closed October29,2026 at02:07 Eastern. Read PLAN.md, REPORT.md and
final-source-evidence.json. The exact317-file candidate has a clean standard
as-cran package check with examples, donttest examples and rebuilt vignettes.
The unchanged exact archive tests ran separately against the check-installed
package with NOT_CRAN=true:5,325passing assertions, zero failures,12existing
warnings and10development-only skips. Focused64 and compatibility263 assertions
are included in that complete count, not additional runs to add to it.

R/spatial_source_report.R is an optional private metadata preflight shared by
likelihood/posterior spatial workflows. It checks supplied projected source CRS,
2D point geometry, linear-unit labels, source/observation IDs and coordinate
columns without transformations, rescaling or re-extraction. It reports nominal
buffer/extent alignment only. Actual valid-cell/kernel-weight coverage and
historical raster-value/extraction-unit lineage remain explicitly unverified.
It does not run inside any fitter or alter public interfaces, prior/numerical
gates, prediction values/attributes or inference eligibility. All existing R,
src, test, DESCRIPTION and NAMESPACE files are unchanged; NEWS alone changes
among old sources, plus the new helper, focused test and generated internal Rd.

unit01's34passes and one missing-raster-CRS error path are retained with original
source/tests/logs. The helper now converts the low-level empty-CRS parse error
into its explicit strict source rejection. No existing gate changed. unit02
passes64assertions with zero failures/errors/warnings/skips; documentation01
passes checkRd with NAMESPACE unchanged; compatibility01 passes263assertions
across spatial targets, preparation, grouped rasters and prediction contracts.

Archive SHA256:
118A69E2A54465CD56EE833A85F80C0B3A379BA0D7AA25084190CC3B40A3FC8A
Workspace:
C:/Users/peterman.73/AppData/Local/Temp/msr-spatial-source-lc238368

Full validation has45resource samples, maximum3actual R interpreters, minimum
13.88GiB available RAM and46.22GiB free disk, zero unknown/foreign processes.
Closeout confirms no R process remains. All five predecessor manifests and
historic pinned evidence pass; no raw chains, tapes or candidate arrays saved.
Read-only verifier from the repository root:

```powershell
& 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' tools/benchmarks/tmb_marginal_likelihood/spatial_source_contract_20261029/verify-candidate.py
```

Never replay this producer or modify this directory after its final manifest.
Later deliberate source edits make the live verifier stale; retain this frozen
coverage and validate changed source separately. Next declared documentation
batch is inference_contract_guide_20261029/PLAN.md and RESUME.md. Its latest
plan refinement reuses the bundled annual Bayesian example as a read-only
consumer, superseding the unexecuted six-point proposal. No fit, chain, new
example dataset or predictive replicate is declared for that guide.

No stage/commit/version bump/merge/push. Authority ends04:45 Eastern; no long
launch after04:00, first wake>=04:15 stops launches, safe checkpoint by04:35 and
delete the one-night heartbeat by04:45. After cutoff checkpoint/delete only.
