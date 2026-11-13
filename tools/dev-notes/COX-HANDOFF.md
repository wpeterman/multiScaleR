# Handoff: Cox prediction metadata under survival 3.8.12 and R 4.6.2

Written 2026-11-12. Read this file first, then `AGENTS.md` in the repository
root.

**Status 2026-11-13: resolved.** The Codex follow-up fixed this
(`R/native_cox_training.R` and changes to `R/conditional_logit.R`,
`R/curve_reference.R`, and `R/fixed_prediction_terms.R`), merged in 0.7.39.
The Cox metadata tests pass with survival 3.8.12. For `clogit()` and `coxph()`
models fitted outside the package with learned transformations, users set
`na.action = stats::na.fail` (see vignettes 07 and 10). The task below is kept
for reference only.

## 1. The task

Eleven failures and errors in five test files occur because multiScaleR can no longer
recover fixed prediction transformations (`poly()`, `scale()`, `ns()` and
nested learners) for `survival::coxph()` and `survival::clogit()` fits. Find
the root cause, fix the package code, and add a focused regression test. Do
not weaken any guard that correctly rejects metadata that truly cannot be
recovered.

These failures are not caused by recent branch work. They occur identically on
commit `bfdc133e3` (0.7.36, branch `inference-dev`) and on `5dd39f5d9` (0.7.37,
branch `docs-restructure`). The package tests last changed for this area on
2026-11-06 (`90ac419df`, release 0.7.17), when they passed.

## 2. Where to work

- Repository: `C:\Users\peterman.73\OneDrive - The Ohio State University\R\Packages\multiScaleR`
  (GitHub `wpeterman/multiScaleR`).
- Create your own worktree and branch from commit `5dd39f5d9` (the tip of
  `docs-restructure`). It holds the newest code, including the related spline
  fix described in section 6, and the test helper scripts. For example:
  `git worktree add <path> -b cox-frame-attributes 5dd39f5d9`.
- Do not edit `.claude/worktrees/docs-restructure` (another session is
  renaming vignettes and editing `_pkgdown.yml` there) or the main checkout
  (branch `inference-dev`). Do not switch their branches.
- The git stash stack is shared with other sessions. Never use bare
  `git stash` or `git stash pop`.
- Your fix touches only `R/` and `tests/testthat/`, so it will merge cleanly
  with the vignette work.

## 3. Rules from the maintainer (Bill Peterman)

- Ask before every commit. Do not push. Do not edit `NEWS` or `DESCRIPTION`.
  Draft the NEWS bullet in your report instead, because the version number is
  assigned when the branches are merged.
- Find the root cause before changing code. Small, reviewable diffs. Follow the
  existing style. Comments say why, not what.
- Every behavior change gets a test in a new file,
  `tests/testthat/test-fix-cox-frame-attributes.R`. Never delete or weaken a
  test to get a pass. Use `skip_if_not_installed()` for suggested packages.
- Roxygen in `R/` uses Rd syntax (`\code{}`), not markdown. Do not edit `man/`
  by hand.
- Prose and comments: American English, active voice, no em-dashes, no
  semicolons.
- R jobs: record the PID of every job you start and only wait on or stop that
  PID. No pattern-based kill commands (`pkill`, `taskkill /IM`,
  `Stop-Process -Name`). A hook blocks them.

## 4. The failing tests

Run with `NOT_CRAN=true`. Messages are identical before and after the 0.7.37
changes.

| File and line | Test | Message |
|:--|:--|:--|
| `test-native-cox-metadata.R:51` | ordinary native Cox rows retain trained transformations and targets | Retained attributes cannot recover the original nested training learner; precompute its predictor or restore its original metadata. |
| `test-native-cox-metadata.R:86` | native Cox checked maps and raw-reference curves keep native scales | same as above (raised inside `predict()` on a raster through `terra::predict`) |
| `test-native-cox-metadata.R:122` | Cox transformation recovery retains unsupported-structure boundaries | Fixed prediction transformations require the retained evaluated training frame. |
| `test-curve-reference.R:210` | retained conditional Cox curves verify raw bases, strata and compound offsets | reference_data does not reproduce the fitted feature: poly(z, 2) |
| `test-fixed-prediction-facets.R:143,144,146` | conditional raw and retained routes share frozen public predictions | `fit` and `se.fit` differ from the oracle (4 failures), then the "Retained attributes" error |
| `test-native-nested-metadata.R:267` | nested conditional choice predictions retain Cox targets and guards | Retained attributes cannot recover the original nested training learner |
| `test-shared-fixed-metadata.R:101` | conditional fixed metadata retains original native expression identities | Retained attributes cannot recover the original nested training learner |

The error sources:

- `R/native_prediction_terms.R:128`, in `.msr_fixed_nested_call()`. It calls
  `makepredictcall()` on the evaluated frame column (`value`) to relearn a
  nested call, so it needs that column's attributes.
- `R/fixed_prediction_terms.R:80`, in `.msr_native_fixed_metadata()`, when the
  retained frame has no `terms` attribute.
- `R/curve_reference.R:127`, where an evaluated reference feature is compared
  with the retained frame column (`.msr_curve_numeric_equivalent()`). Check
  whether that comparison includes attributes.

Callers that pass a Cox frame: `.msr_native_fixed_metadata()` uses
`model$model` for `coxph`. `.msr_clogit_contract()` in `R/conditional_logit.R`
calls `.msr_fixed_prediction_terms(tt, m$model)`.

## 5. What is known about the cause

Probe (survival 3.8.12, R 4.6.2, script below):

```r
library(survival)
set.seed(1)
d <- data.frame(time = rexp(60), status = rbinom(60, 1, .7), z = rnorm(60), x = rnorm(60))
f <- Surv(time, status) ~ stats::poly(z, 2) + x
m <- coxph(f, data = d, model = TRUE)
names(attributes(m$model[[2]]))                       # dim, dimnames
names(attributes(model.frame(m)[[2]]))                # dim, dimnames
mf <- model.frame(f, data = d)
names(attributes(mf[[2]]))                            # dim, dimnames, coefs, degree, class
names(attributes(na.omit(mf)[[2]]))                   # dim, dimnames
m2 <- coxph(f, data = d, model = TRUE, na.action = na.pass)
names(attributes(m2$model[[2]]))                      # dim, dimnames, coefs, degree, class
attr(terms(m), "predvars")                            # stats::poly(z, 2, coefs = list(...)) is learned
```

So:

1. The Cox model frame (`m$model`, and `model.frame(m)`) has lost the class
   and recipe attributes of matrix columns (`poly`, `ns`, `bs`, and the
   attributes of `scale()`). It keeps only `dim` and `dimnames`.
2. Running a model frame through `na.omit()` strips them the same way, even
   when no row is removed, because `[.data.frame` row subsetting drops them.
   `coxph()` with `na.action = na.pass` keeps them. `model.frame()` alone keeps
   them.
3. The saved prediction call, `attr(terms(m), "predvars")`, still holds the
   learned recipe for a top-level `stats::poly()`. Nested learners may not be
   learned there, which is why `.msr_fixed_nested_call()` reads the frame.

Not yet known: which environment change started this. Two candidates arrived
together: R 4.6.2 (released 2026-10-29) and survival 3.8.12. The survival
NEWS for 3.8-10 says "Fix github issue 341 (model.frame.coxph)". Settle this
before choosing a fix, for example by installing an older survival (3.8-3 or
earlier) into a private library (`install.packages(..., lib = <scratch>)`,
then put it first in `.libPaths()`) and rerunning the probe and one failing
test. Do not change the user's library. GitHub Actions history
(`gh run list --repo wpeterman/multiScaleR`) may also show when these tests
started failing.

## 6. A related fix already made (precedent)

Commit `5dd39f5d9` fixed the same symptom for fixed-knot splines with
`backend = "conditional_logit"`. See `.msr_spline_predvars_match()` in
`R/response_spline.R` and its use in `.msr_compile_spline_design()`. When a
frame column is a bare matrix with no basis attributes, it compares the fixed
recipe with the saved `predvars` call. It accepts only literal constants and
never evaluates caller code. A column that still carries attributes must match
them directly. The test is `tests/testthat/test-fix-clogit-splines.R`, with a
step-selection fixture in `tests/testthat/helper-step-selection.R`.

## 7. Approaches to evaluate (not prescribed)

- Use the learned call in `predvars` when it is already complete, so the frame
  attributes are not needed. Check `.msr_fixed_normalize_call()` and
  `.msr_needs_nested_metadata()` to see why that path is not taken for these
  fits.
- For nested learners, rebuild an attribute-preserving frame for metadata
  recovery only. One route is to re-evaluate the fitted terms on the recovered
  training data with `na.action = na.pass`, only when the fit dropped no rows
  (`is.null(m$na.action)`) and the row names match the retained frame exactly.
  Never use it to change fitted values.
- Restore the missing attributes from the learned `predvars` call onto a copy
  of the frame column, again only with literal constants.

Whatever you choose, the guards that refuse two learners in one sum,
unsupported structures, penalized or time-varying Cox models, and genuinely
missing metadata must keep refusing them. `test-native-cox-metadata.R:122`
("retains unsupported-structure boundaries") is a good check that guards still
fire.

## 8. Running code and tests

- R: `C:/Program Files/R/R-4.6.2/bin/Rscript.exe`. Use Git Bash or PowerShell.
  `Rscript -e` with pipes breaks under Windows quoting, so write script files.
- The user library holds an old multiScaleR 0.7.2. Always use
  `devtools::load_all()` from your worktree, never `library(multiScaleR)`.
- `tools/dev-notes/scripts/run_tests.R` loads the source and prints one summary
  row per test file:
  `"$RS" --vanilla tools/dev-notes/scripts/run_tests.R test-native-cox-metadata.R ...`
  (run from the worktree root, with `NOT_CRAN=true`).
- The compiled code builds on first `load_all()` (about 2 to 3 minutes). To
  skip that, copy `src/*.o` and `src/multiScaleR.dll` from the
  `docs-restructure` worktree into your `src/` (the C++ sources are identical)
  and `touch` them.
- Ignore the "glmmTMB was built with TMB 1.9.21" warning. It is a local build
  mismatch.
- Tests that start PSOCK workers load the installed package, so under
  `load_all()` they fail against the old 0.7.2 in the user library. Judge
  those under `R CMD check` only.

Validation before you report:

1. The five files above pass.
2. No regressions in the files that exercise the same code:
   `test-conditional-logit.R`, `test-amt-conditional-logit.R`,
   `test-kernel-splines.R`, `test-fix-clogit-splines.R`,
   `test-cv-choice-spatial.R`, `test-native-raw-transform.R`,
   `test-fixed-kernel-radius.R`,
   `test-unmarked-trained-prediction.R`, `test-gdr-trained-prediction.R`.
   Grep `tests/testthat` for `coxph`, `clogit`, `poly(` and `scale(` to find
   others.
3. If feasible, the full `devtools::test()` (about 50 minutes, run it in the
   background and record the PID).

## 9. Report back

- The root cause (which package or R change, with evidence), the file and line
  of each fix, and why it is correct.
- Files changed and the new test file.
- Any user-visible change in numbers, messages or defaults, and a drafted NEWS
  bullet.
- Test summary rows (file, tests, failed, errors, skipped).
- Do not commit until Bill says so.
