# Phase63: opt-in matching of raster and prepared neighborhoods

## FINAL checkpoint, October 28 at 14:10 Eastern

Today's bounded batch is complete ahead of the 15:30 logoff. Automation
`multiscaler-spatial-workflows-until-3-30pm` has been DELETED successfully. No
task-owned worker remains active, and no further launch is scheduled. Other R
processes run the owner's `latent_gaussian_20261028/pilot-null.R`, started after
both projection timings completed; leave them untouched. Final process evidence
is `final-process-checkpoint.json`. Historical progress notes below are retained.

Current authoritative evidence is REPORT.md. Combined attempt2 is complete:
0 errors,0 warnings,0 notes, including rebuilt vignettes and 4,763 installed-test
assertions. Twelve existing fixture warnings and ten development-only skips are
recorded separately. The final source audit matches all305 checked source files,
with no additional package sources. New opt-ins include_fixed/match_buffer remain
FALSE by default; all safeguards/starts/Bayesian edits are preserved. No commit.

Both performance pairs are accepted and frozen. The saved N-mixture repeated
snapshot workload is3.21s live versus0.03s cached plus0.19s preparation, peak
RSS1,246,011,392 ->870,002,688bytes; cold startup was slower. Equivalent mean/SQ
transforms over the existing2,500-cell grid are9.95s point preparation versus
0.79s matched raster, peak RSS951,840,768 ->789,594,112bytes and allocations
524,439,128 ->8,260,024bytes. All-cell differences are<=3.89e-16; both contention
records are empty. These are separate single-pair scopes, not large-map claims.
Rerunning either accepted pair or the combined check is unnecessary absent a
new source change, concrete failure or new validation question.

For a later authorized turn, first inspect processes and rerun
`audit-checked-source.py 2` using the Python3.12 path below. Preserve owner drift
and qualify any now-stale combined coverage. Next priorities are explicit
source/CRS/unit and coverage diagnostics, a representative larger existing-map
performance gate, and a separately defined statistical calibration design.
Do not silently add scale uncertainty, reconstruct opaque native frames, change
metric definitions, or treat relative-selection maps as choice probabilities.

## Historical working handoffs

Continue only through15:30 America/New_York October28,2026. No long launch after
14:45; checkpoint by15:10 and stop/delete today's automation before15:30.
Automation remains multiscaler-spatial-workflows-until-3-30pm; its first handoff
is phase62/RESUME.md. Preserve all owner/Bayesian work, defaults, independent
starts/cold audits, explicit failures and uncertainty safeguards. No commit,
stage, version bump, merge or push. No new datasets/optimizer panels.

Phase62 accepted final-focused-3:351 assertions,0fail/error/warning/skip. Verified
fixed-radius native point/snapshot and unmarked/GLM/coxph raster preparation,
native predictor discovery, training-row identity and spatial-only clamping.
Native relative-risk map contract is tested and documented. Defaults unchanged.

Next bounded feature: match_buffer=FALSE by default in kernel_scale.raster().
Opt-in TRUE uses the retained preparation's max_D polygon-intersection support for kernel
means and weighted surface metrics, with finite-cell normalization. Retain
center-cell NA masks so the raster domain is not filled by neighbors. Reject fixed
radii above max_D because capping a center-distance circle does not reproduce the
actual truncated polygon-intersection footprint. Initial scope FFT/na.rmTRUE
on square projected grids; reject incompatible requests explicitly. Current
default pct_wt window remains numerically unchanged. Binned fitting remains an
approximation; do not advertise strict binned/exact equality or uncertainty
calibration. Weighted SQ stays population weighted spread; hard-radius SQ stays
sample SD. Validate on established surface/landscape fixtures and existing cell
centers, including boundaries/missing masks. No optimizer work is justified.

Unrelated Bayesian workers19420 and70564 are active as of12:59; RStudio recently
closed. Do not benchmark or full-check concurrently. Source-bound non-timing
validation is allowed under conservative resource gates. Record accepted gains,
allocation-only evidence and pending time/peak/package-check gates separately.

13:20handoff: first strict projection tests failed because the simulated fixture
had no projected CRS label; the test copy now explicitly labels its known meter
grid. Second tests exposed boundary cells: kernel_prep retains every polygon-
intersecting cell, not only centers within max_D. Matching weight support now
uses the same sf::st_buffer/exact_extract construction on a centered template.
No fitting behavior or tolerance was changed. Both failed logs are retained.
Focused run3 is active via phase63/validate.R 3; inspect focused-3.log/process.
Package and unrelated terradish checks are active; do not duplicate them.

The pending phase61 time/peak pair completed successfully in a brief uncontended
window. Both contention logs empty; strict source/fit/native gates pass. Across
200 predictions3.21s ->0.03s plus0.19s preparation; sampled process peak1.246GB ->
0.870GB. Cold process wall14.687 ->16.078s was slower. Read the appended phase61
REPORT section and spatial-audit.rds. Do not rerun accepted outputs. Source
matching is a later opt-in and does not invalidate those frozen default-path
measurements. Report timings/peak and cold startup separately.

## 13:55 current handoff and exact commands

Consolidated current evidence: REPORT.md. Do not use the earlier process
inventory as current ownership proof. Fixed-prediction final validation is 351
assertions, matching focused run 3 is 75, and broad focused validation is 1,454
with three retained existing binomial fixture warnings, zero failures/errors.
The missing new-export inventory entry was corrected; focused inventory test is
21 assertions, zero failures/errors/warnings. The first known-stale full attempt
was stopped safely and retained. Fresh full attempt 2 is active. No duplicate
checks or accepted benchmark reruns; no package source edits while it runs.

Inspect progress in PowerShell from the repository root:

```powershell
Get-Date -Format o
Get-Content tools/benchmarks/joint_scale_assessment/phase63/check-controller-2.json
Get-Content tools/benchmarks/joint_scale_assessment/phase63/check-controller-run-2.log -Tail 30
```

After completion inspect `check-release-2/check-summary.txt`,
`check-release-2/R-CMD-check.log` and `check-release-2/installed-tests.log`.
Run `audit-checked-source.py 2` with the Python 3.12 executable to audit live
source identity and additional source files. Only advertise exact combined
coverage if the manifest and new-file scan match. The controller's state/logs
retain root and descendant process creation times; use those for ownership if
an unfinished task must stop. Controller itself has a 15:10 Eastern cutoff and
resource/runtime reserves. Preserve all owner jobs.

At the first run at or after 15:10, finish REPORT.md and this handoff, stop new
launches, preserve completed atomic outputs and delete automation
`multiscaler-spatial-workflows-until-3-30pm` before 15:30. Do not create another
nightly automation or recur after this date. No staging, commits or version bump.

After a successful combined check and before 14:45, one bounded serial matched-
projection pair is prepared in `projection-pair.py`. It recreates the exact
existing 50x50, resolution10, seed99 raster fixture, without its unrelated model
fit, and evaluates unbinned Gaussian mean/SQ at all 2,500 cell centers, with the
same 150-unit buffer and independent scales60/100. Reference uses point
preparation; candidate uses the matching raster branch. This is a numerical
workflow comparison, not a new scientific dataset or optimizer panel. The script
refuses other R surfaces, applies RAM/disk/cutoff gates and stops only verified
owned processes. Run once using Python3.12, then audit with R4.6.1:

```powershell
& 'C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe' tools/benchmarks/joint_scale_assessment/phase63/projection-pair.py
& 'C:/Program Files/R/R-4.6.1/bin/x64/Rscript.exe' --vanilla tools/benchmarks/joint_scale_assessment/phase63/audit-projection.R
```

Require empty competing-process records, successful output equivalence at1e-10
and matching source hashes before accepting timings. Report native allocations
and sampled process peaks separately from Rprofmem; retained reference includes
reusable point buffers while raster output has a different representation.
Do not overwrite partial or accepted pair directories or rerun a failed mode
without inspecting the exact failure first.

14:05 update: attempt 2 finished **0 errors, 0 warnings, 0 notes**, vignette-
enabled, 16m31.3s. Installed tests:4,763 pass,0fail,12 retained fixture warnings,
10 development-only skips. The post-check manifest audit is exact (305 files,
no additional package sources). The same-grid projection pair is now active
via projection-pair.py, output projection-pair.log; inspect before duplicating
anything. Its per-mode atomic results and ownership records are under
projection-reference-1/ and projection-candidate-1/. No package source edits
since the combined check freeze.
