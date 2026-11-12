# Handoff: multiScaleR documentation and code update

Written 2026-11-12. Read this file first, then `AGENTS.md` in the repository
root, then `vignette-relocation-map.md` in this folder.

## 1. Where the work is

- Worktree: `<repo>/.claude/worktrees/docs-restructure`, branch
  `docs-restructure`, cut from `inference-dev` at `bfdc133e3` (version 0.7.36).
  The main checkout stays on `inference-dev`; do not switch its branch.
- **Nothing is committed.** Every change below is in the working tree.
- The git stash stack is shared with other sessions. Never use bare
  `git stash` or `git stash pop`.
- Reusable scripts are in `tools/dev-notes/scripts/` (this folder is excluded
  from the built package by `.Rbuildignore`). They hard-code a scratch
  directory in a `scr <-` line near the top. Point that at your own scratch
  directory before use.
- R is `C:/Program Files/R/R-4.6.2/bin/Rscript.exe`. The user library holds an
  old multiScaleR 0.7.2. To render vignettes against the branch, install the
  worktree into a private library and put it first in `.libPaths()`:
  `R CMD INSTALL --no-multiarch --no-docs --library=<scratch>/rlib .`
  **Reinstall after the code fixes below, before re-rendering vignettes.**
- For tests, `scripts/run_tests.R` uses `devtools::load_all()` and prints one
  summary row per test file.

## 2. Standing rules from the maintainer (Bill Peterman)

- **Ask before every commit.** Bump the version by 0.0.1 per commit and add a
  `NEWS` entry under the new 0.7.x number (plain `NEWS`; it becomes `NEWS.md`
  only at the 1.0.0 release, as a rename-only commit).
- Run `devtools::document()`, `devtools::test()`, and
  `devtools::check(args = "--no-manual")` before a commit.
- Never overwrite his edits. If he hands a file back, diff before editing.
- R jobs: record the PID of every job you start and only wait on or stop that
  PID. No pattern-based kill commands.
- Writing: no em-dashes, no semicolons in prose, colons sparingly, active
  voice, American English, "scenario" for a simulation design combination.
  The `peterman-research-suite` skills are the standard (academic-voice-style,
  r-doc-style, r-vignette-writing, humanizer as the last pass).
- `terra` first; `sf` only where required. Do not invent references.
- Do not edit `.Rd` files by hand. Add tests for code changes. Never weaken a
  test to get a pass.
- Keep his name, email, and project details out of requests to outside
  services.

## 3. State of the vignettes

**Bill reviewed the rendered vignettes and approved them on 2026-11-12.**
Fifteen vignettes exist: 13 rewritten, 2 new (`step-selection.Rmd`,
`spatial-response-fields.Rmd`). All 15 built cleanly on 2026-11-11 (about four
minutes in total). Rendered copies are in `.claude/pilot-review/` (ignored by
git). Two lines changed after his review: the Guide's `plot_kernel()` call now
uses `shape = 5`, and the projection vignette states the new `pct_mx` bounds.

Open content choices he was told about and did not ask to change (leave as
is unless he says otherwise): the projection vignette dropped the
biased-sampling simulation; several vignettes exceed their line targets; the
spatial-fields and distance-removal examples use seeds chosen after trying
several; three Bayesian remedies are standard Stan practice, not package text.

## 4. Code changes: what is done

Targeted tests passed for each item when it was made. No full `devtools::test()`
or `R CMD check` has been run on the combined tree.

| Change | Files | Tests |
|:--|:--|:--|
| `plot_marginal_effects(add_label = FALSE)`: captions are opt-in, including the Bayesian caption | `R/marginal_effects.R`, `R/bayes_marginal_effects.R` | 9 existing files updated to request the caption; new default test in `test-native-marginal-effects.R` |
| `pct_mx` is range-based: limits are `min - pct_mx * range` and `max + pct_mx * range`; valid range is greater than -0.5 up to 0.99 | `R/scale_center_raster.R`, `R/kernel_scale.raster.R` | `test-pct-mx-range.R` |
| `shape` replaces `beta` in `kernel_dist()` and `plot_kernel()`; `beta` warns for one release; both with different values is an error | `R/kernel_dist.R`, `R/plot_kernel.R`, callers in `R/scale_vars.R`, `R/surface_metrics.R`, `R/kernel_scale.raster.R` | `test-shape-argument.R` (see section 5, item 1) |
| `select_cv()` scores on the metrics shared by all candidates; result gains `selection$metrics` and `selection$metrics_dropped` | `R/cv_comparison.R`, `R/cv_distribution_scores.R` | `test-fix-select-cv-mixed.R` |
| `landscape_var()`: warns when `base` is given with `msidi` or `msiei` (they use the natural log); `rpr` without `classes_max` is an immediate error; help corrected | `R/scale_vars.R` | `test-fix-landscape-base-rpr.R` |
| Un-exported `audit_gdr_integration()` and `score_simulation()` (`@noRd`); `vcov.multiScaleR` is an S3 registration only | `R/audit_gdr_integration.R`, `R/score_simulation.R`, `R/profile_glm_uncertainty.R`, `R/prepare_gdr_histories.R`, `_pkgdown.yml` | none yet |
| `sf` to Imports, `pscl` to Suggests, `graphics` and `grDevices` declared | `DESCRIPTION` | none yet |
| Absolute build paths in `inst/extdata/bayesian-annual-example.rds` replaced with `R/...` and `src/...`; rest of the object verified identical | that file | `test-fix-shipped-data-paths.R` (written, not yet run) |
| Land-cover breaks in `data-raw/generate_vignette_cache.R` moved inside the 0 to 1 range of `land2` | that file | cache not regenerated |

`gdr_start_values()` **stays exported.** Bill confirmed this on 2026-11-12
(it replaces his earlier decision to un-export it), because
`?fit_gdr_response_field` and `?make_gdr_refitter` tell users to call it.

## 5. Code changes: what is left

**Status 2026-11-12 (second session): items 1 to 11 are done, uncommitted.**
Each has a test (`test-fix-clogit-splines.R`, `test-fix-optimizer-start.R`,
`test-fix-field-geometry.R`, `test-fix-profile-backends.R`, additions to
`test-fix-landscape-base-rpr.R`). Fix status per issue is in the "Code issues"
table of `vignette-relocation-map.md`. `devtools::document()` was run (NAMESPACE
diff as expected), `NEWS` has a 0.7.37 entry, and `DESCRIPTION` is 0.7.37.
Full `devtools::test()` and `devtools::check(args = "--no-manual")` were run.
The remaining test errors (Cox prediction metadata in `test-native-cox-metadata.R`,
`test-curve-reference.R`, `test-fixed-prediction-facets.R`,
`test-native-nested-metadata.R`, `test-shared-fixed-metadata.R`) occur
identically at `bfdc133e3` and come from survival 3.8.12 stripping basis
attributes in `coxph()` frames. The two RTMB-dependent checks now skip when
RTMB is absent. Check NOTE: a hidden `.git` file (worktree artifact). Built
tarball: 8.3 MB, mostly `bayesian-annual-example.rds` (2.9 MB, unchanged) and
the other example and cache files. The fixed-kernel default search is
unchanged (warning only), pending Bill's decision on whether `kernel = "fixed"`
should default to screened starts. The original item list follows for
reference.

Evidence and reproduction recipes for every issue: the "Code issues found
during the rollout" table in `vignette-relocation-map.md`, and
`scripts/verify_issues.R` with its output `scripts/verify_issues.log`.
`scripts/fix_brief.md` holds the rules given to the fix agents.

1. **Known failing test.** `tests/testthat/test-shape-argument.R` has a test
   named "package-internal callers that still pass beta do not warn users".
   It tested a guard that has been removed (internal callers now use `shape`).
   Delete that test, then rerun `test-shape-argument.R`, `test-kernel-math.R`,
   `test-planar-kernel-reporting.R`, `test-surface-metrics.R`,
   `test-landscape-metrics.R`, and `test-coverage-expansion.R` (the last needs
   `NOT_CRAN=true`).
2. **Unreviewed partial edit.** `R/response_spline.R` has 48 added lines from
   an agent that was cut off while working on fixed-knot splines with
   `backend = "conditional_logit"`. It has no test and has not been read or
   run. Read the diff (`git diff -- R/response_spline.R`). Keep it only if
   you can show it is correct; otherwise restore the file with
   `git checkout -- R/response_spline.R` after confirming nothing else there
   is wanted.
3. **Fixed-knot splines with `conditional_logit`.**
   `clogit(case_ ~ ns(forest, knots = 0, Boundary.knots = c(-3, 3)) + sl_km + strata(step_id), model = TRUE)`
   then `multiScale_optim(backend = "conditional_logit")` stops with "Native
   dynamic spline frame does not match its fixed basis recipe." The help
   (`R/multiScale_optim.R`, about lines 711 to 727) says this is supported,
   and a passing test in `test-kernel-splines.R` ("fixed splines retain native
   mixed count and choice-set targets") builds its model differently. Find
   why the user-level route fails. Fix the code, or correct the help. Add a
   user-level test. New test file: `test-fix-clogit-splines.R`.
4. **Optimizer returns its start values silently.** Two cases, both end at
   415 m for every layer with `NA` standard errors and no package warning:
   (a) `kernel_prep(kernel = "fixed")` with the default start (logLik -175.04;
   `start_strategy = "screen"` reaches 410 and 724 m, logLik -154.63);
   (b) a glmmTMB field model fitted from glmmTMB's default starts
   (non-positive-definite Hessian). Find the cause in each case. Make the
   fixed-kernel default search work if that fits how hard radii are handled
   elsewhere. In all cases warn once, name the remedy, and record it in the
   fit's diagnostics so `summary()` shows it. Do not change results for fits
   that already converge. New test file: `test-fix-optimizer-start.R`.
   If you change fixed-kernel behavior, the Guide's "Choose a start strategy"
   section uses this stall as its example and needs its prose revisited.
5. **`response_field_geometry()` rejects real fits.**
   `R/response_field_geometry.R:117` requires `object$backend == "refit"`, and
   a default-backend fit has no `backend` element. Existing tests pass only
   because fixtures set it by hand. Fix, and test with a real (small) field
   fit built through `multiScale_optim()`. glmmTMB needs explicit starts on
   meter coordinates, for example `start = list(theta = c(log(0.5), log(500)))`.
   New test file: `test-fix-field-geometry.R`. Once fixed, the spatial-fields
   vignette can use the geometry route for its map in place of the grid of
   points (optional).
6. **Unclear error at new field positions.** `predict()` on a field fit at
   new positions without `allow.new.levels = TRUE` gives "subscript out of
   bounds". Give a clear message naming `allow.new.levels = TRUE` and
   `re.form = NA`.
7. **`profile_sigma()` defaults fail on a `profile_glm` fit** ("sigma values
   must remain within the fitted bounds"). Clip the default grid inside the
   bounds; keep an informative error for a user-supplied range. Check
   `summary(fit, profile = TRUE)` too. New test file:
   `test-fix-profile-backends.R`.
8. **Profiling a `conditional_logit` fit** errors unhelpfully, after a stray
   "first element used of 'length.out'" warning. Support profiling if it can
   be checked against a grid of refits; otherwise stop early with one clear
   message. Fix the stray warning either way.
9. **Help wording in `R/multiScale_optim.R`.** Say that `profile_glm` needs a
   `glm()` fit for Gaussian models (`lm()` is rejected with a clear error).
   Remove the line `#' @importFrom pscl zeroinfl` (about line 1053), since
   `pscl` is now in Suggests and `zeroinfl` is never called. Confirm tests
   that use `pscl::zeroinfl` (`helper-fixtures.R`, `test-internals-and-helpers.R`)
   skip cleanly when `pscl` is absent.
10. **Not rerun.** An agent reported that `multiScale_optim()` on a
    distance-removal field fit returned a sigma standard error of 3.4e-30 and
    that its map step failed. The run takes about seven minutes. Reproduce
    before deciding anything.
11. **Rough edges noted, lower priority:** `amt`-based and `RTMB`-based routes
    were not exercised because `RTMB` and `tmbstan` are not installed here;
    other diversity metrics also accept a `base` they ignore without warning.

Then, for the code-change commit:

- `devtools::document()`. Roxygen markdown is still off for this commit.
  Review the `NAMESPACE` diff: `export(audit_gdr_integration)`,
  `export(score_simulation)`, `export(vcov.multiScaleR)`, and
  `importFrom(pscl,zeroinfl)` should disappear; `S3method(stats::vcov,multiScaleR)`
  and the print method for `msr_simulation_score` must remain.
- Full `devtools::test()` and `devtools::check(args = "--no-manual")`.
- `NEWS` entry under 0.7.37 and `DESCRIPTION` version 0.7.37. User-visible
  items: marginal-effect captions off by default; `pct_mx` rule and bounds
  (changes projected values when `pct_mx` is not 0); `shape` with `beta`
  deprecated; `select_cv()` on mixed metric sets; `msidi`/`msiei` `base`
  warning; `rpr` requires `classes_max`; two functions no longer exported;
  `sf` no longer attached with the package (user scripts that relied on that
  need `library(sf)`); plus whatever items 3 to 8 change.
- Show Bill the diff summary and ask before committing.

## 6. Adding the vignettes to the update and the website

Bill's decisions (all made explicitly; details in the project memory file
`docs-redteam-vignette-plan.md`): number the vignettes by file-name prefix and
index title; quick start first, then the Guide, then by complexity; release
site only, deployed from `master`.

Do these after the code-change commit, as one or two commits, asking first.

1. **Reinstall and re-render** all 15 vignettes against the fixed code
   (`scripts/render_all.R`). Check that numbers quoted in prose still match
   output. Known places that move: the `pct_mx` table in
   `spatial_projection_clamping.Rmd` (rows for -0.1 and 0.1 now follow the
   range rule, which the prose already describes), and anything touched by
   items 3 to 8 above. Restrictions tables to update if behavior changes:
   `step-selection.Rmd` (profile intervals, splines),
   `feature-compatibility.Rmd` (spline row for `conditional_logit`,
   `profile_sigma()`), `spatial-response-fields.Rmd` (geometry route,
   starting values), `cross-validation.Rmd` (it works around the `select_cv()`
   metric mismatch, which is now fixed).
2. **Consolidate the caches.** Five scripts `data-raw/cache_*.R` write five
   files `inst/extdata/*-cache.rds` (about 0.7 MB together), next to the
   older `vignette_cache.rds` (0.7 MB). Bill wants one regeneration script.
   Fold them into `data-raw/generate_vignette_cache.R` (one script may still
   write several small files). The `landscape` and `surface` entries of
   `vignette_cache.rds` are no longer read by any vignette and can go; the
   `quickstart` entry is still used, and the Quick Start relabels its stored
   call in a hidden chunk, which can be removed once the cache is rebuilt with
   matching object names. Check that no vignette reads `inst/extdata/*.RData`
   (`.Rbuildignore` excludes them, so they do not ship) and report the built
   tarball size against CRAN's 5 MB guidance.
3. **Rename and number.** File names must start with a letter. Bill approved this
   exact list of file names and the numbered index titles on 2026-11-12:

   | New file | Current file | Index title |
   |:--|:--|:--|
   | `v01-quickstart` | `quickstart` | 01. Quick start |
   | `v02-user-guide` | `multiScaleR_Guide` | 02. User guide |
   | `v03-map-predictions` | `spatial_projection_clamping` | 03. (current title) |
   | `v04-landscape-metrics` | `landscape_metric_covariates` | 04. (current title) |
   | `v05-surface-metrics` | `surface_metric_covariates` | 05. (current title) |
   | `v06-year-specific-rasters` | `year-specific-rasters` | 06. (current title) |
   | `v07-fitting-methods` | `feature-compatibility` | 07. (current title) |
   | `v08-cross-validation` | `cross-validation` | 08. (current title) |
   | `v09-estimates-and-intervals` | `inference-contracts` | 09. (current title) |
   | `v10-step-selection` | `step-selection` | 10. (current title) |
   | `v11-spatial-fields` | `spatial-response-fields` | 11. (current title) |
   | `v12-distance-removal` | `distance-removal-audits` | 12. (current title) |
   | `v13-bayesian-scales` | `bayesian-scales` | 13. (current title) |
   | `v14-bayesian-prediction` | `bayesian-prediction` | 14. (current title) |
   | `v15-bayesian-computation` | `bayesian-computation` | 15. (current title) |

   Each index title is the zero-padded number followed by the vignette's
   current title (for 01 and 02, "Quick start" and "User guide").

   Use `git mv` so history follows. Zero-pad the numbers in
   `\VignetteIndexEntry` (R sorts `browseVignettes()` by title; CRAN and the
   installed index sort by file name). Check whether `html_vignette` warns
   when the YAML `title` differs from the index entry and keep them equal if
   it does. Then update every `vignette("...", package = "multiScaleR")`
   cross-link in all vignettes (each has a "See also" list; the Quick Start
   has a "Where to go next" table), every `\code{vignette(...)}` in `R/`
   roxygen, the README, and the cache file names if they embed vignette names.
   Grep for each old name afterward.
4. **Website (`_pkgdown.yml`).** Update `articles:` to the new file names and
   add the two new articles; pkgdown's check fails unless every vignette is
   listed. Suggested groups: Getting started (01, 02); Covariates and maps
   (03 to 06); Fitting and evaluation (07 to 09); Specialized models (10 to
   12); Bayesian scale estimation (13 to 15). Update the navbar "Get started"
   link from `articles/quickstart.html` to the new file. Under `reference:`,
   the two un-exported functions are already removed; add any new help topics
   from section 7. Run `pkgdown::check_pkgdown()` and build locally with
   `pkgdown::build_site()` to confirm article order and that every
   cross-link resolves. pkgdown renders only `NEWS.md`, so the changelog page
   stays absent until 1.0.0. Do not deploy or push; the site publishes from
   `master` when Bill merges.
5. **README.** It cites 4 of the old 13 vignettes. Point it at the numbered
   set.
6. **`NEWS`.** Add one user-facing "Changes since 0.7.4" summary block
   (Bill's decision), covering the new vignettes and structure.
7. Delete `tools/dev-notes/vignette-relocation-map.md` only when every row
   reads done, including the help-page rows.

## 7. Help-page phase (after the vignettes land)

Decided by Bill; none started:

- Enable roxygen markdown globally (`Roxygen: list(markdown = TRUE)`),
  regenerate, and review the full `man/` diff. Fix unescaped `%` in
  `R/multiScale_bayes.R` (about lines 90 and 376).
- `@noRd` for the 25 internal help pages.
- Restructure the long pages: `multiScale_optim`, `multiScale_bayes`,
  `predict.multiScaleR`, `cross_validate`, `kernel_scale.raster`.
- Move the content listed in `vignette-relocation-map.md` (table "Rollout of
  the remaining vignettes" and the pending rows above it) onto the named help
  pages.
- Fix the profile-interval text: `summary(fit, profile = TRUE)` is the route.
  Complete the `\value` of `summary.multiScaleR` (5 of 15 elements named).
- A to Z help index: R builds the letter bar only above 100 indexed topics.
  Add aliases for undocumented print methods, one page per returned object,
  and topic pages for kernels, fitting methods, prediction quantities, and a
  glossary (`?multiScaleR-glossary`, to be linked from every vignette).
  Verify the bar in the built `00Index.html`.
- A runnable example for every export. Rewrite `?multiScaleR` and the
  `DESCRIPTION` text.
- Test timing pass with `skip_on_cran()`; no test removed or weakened.

## 8. Files in this folder

- `HANDOFF.md`: this file.
- `vignette-relocation-map.md`: where removed vignette content goes; code
  issue verdicts with evidence; cache consolidation note.
- `validation-status.md`: development-status text moved out of the vignettes,
  verbatim, grouped by source file.
- `scripts/`: `verify_issues.R` and its log, `run_tests.R`, `render_one.R`,
  `render_all.R`, `html_inspect.R`, `prose_only.R`, `check_cv_spline.R`, and
  the two agent briefs.
