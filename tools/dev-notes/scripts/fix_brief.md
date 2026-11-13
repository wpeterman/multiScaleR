# Code-fix brief (multiScaleR, branch docs-restructure)

You are fixing confirmed defects in the R package `multiScaleR` for its
maintainer. Several agents work in the same worktree at once on DIFFERENT
files. Stay inside the files your task assigns.

## Where things are

- Worktree (work ONLY here; never `cd` to the parent repository; never switch
  branches): `C:\Users\peterman.73\OneDrive - The Ohio State University\R\Packages\multiScaleR\.claude\worktrees\docs-restructure`
- Scratchpad (SCR): `C:\Users\peterman.73\AppData\Local\Temp\claude\9e3b7d3898bd820b\459cb932-8333-43df-b0da-7d9a2422aa59\scratchpad`
- Evidence for every issue: `tools/dev-notes/vignette-relocation-map.md`,
  section "Code issues found during the rollout", and the reproduction script
  `SCR/verify_issues.R` with its output `SCR/verify_issues.log`. Read the
  parts for your issues first and reuse their recipes.
- Project rules: `AGENTS.md` in the worktree root. Read it.

## Rules

- Find the root cause before changing code. Read the surrounding code and the
  existing tests. A fix that only hides a symptom is not acceptable. If the
  right behavior is a scientific or statistical judgment you cannot settle
  from the code and documentation, do NOT guess: make the smallest safe change
  (for example a clear error or warning), and report the open question.
- Small, reviewable diffs. Preserve public interfaces unless your task says
  otherwise. Do not silently change numerical results. If a fix changes a
  number a user would see, say exactly which and why in your report.
- Follow the existing code style of each file. Add moderate comments that say
  why, not what.
- Every behavior change gets a focused test in a NEW file
  `tests/testthat/test-<topic>.R` named in your task (do not append to
  existing test files unless an existing expectation must change because the
  contract changed; say so). Tests must fail before your fix and pass after.
  Keep them fast (seconds). Use `skip_if_not_installed()` for suggested
  packages. Never delete or weaken an existing test to get a pass.
- Roxygen comments in `R/` are the documentation source. Update them when
  behavior or wording changes. This package's roxygen uses Rd syntax
  (`\code{}`, `\link{}`), not markdown. Do NOT run `devtools::document()`, do
  NOT edit `man/`, `NAMESPACE`, `DESCRIPTION`, `NEWS`, `_pkgdown.yml`,
  `vignettes/`, `src/`, or `inst/`. The lead regenerates documentation once.
- Documentation text: American English, no em-dashes, active voice, concise.
  Do not invent references or claims.
- No git write commands (no commit, stash, checkout, reset, add).
- Do not install packages or modify `SCR/rlib`.
- Job discipline: launch each R job so you know its PID
  (`cmd > log 2>&1 & echo "PID $!"; wait $!`), and only wait on or stop that
  exact PID. Never use `pkill`, `killall`, `taskkill /IM`, `Stop-Process -Name`.
- If you need to touch a file outside your assignment, stop and report it
  instead.

## Running code and tests (Git Bash)

```bash
cd "<worktree>"
SCR="C:/Users/peterman.73/AppData/Local/Temp/claude/9e3b7d3898bd820b/459cb932-8333-43df-b0da-7d9a2422aa59/scratchpad"
RS="/c/Program Files/R/R-4.6.2/bin/Rscript.exe"
# run test files against the worktree source (uses devtools::load_all), one summary row per file
"$RS" --vanilla "$SCR/run_tests.R" test-a.R test-b.R
```

For exploration, write a script `SCR/fix_<topic>_explore.R` that starts with
`suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))` so it
uses the worktree source (the user library holds an old version). `Rscript -e`
with pipes breaks under Windows quoting, so use script files. The package's
compiled code is already built; do not trigger a rebuild by touching `src/`.
Other agents run `load_all()` at the same time, so if a load fails oddly once,
retry before diagnosing.

Before you finish, run your new tests plus every existing test file that
covers the code you changed (find them with grep), and report the summary
rows.

## Report back (concise)

- For each issue: root cause (file and line), what you changed, and why that
  is the right behavior.
- Files edited and new test files.
- Any user-visible change in numbers, messages, or defaults (for the NEWS
  entry).
- Test summary rows (file, tests, failed, errors, skipped).
- Help-text changes needed in files you were not allowed to edit.
- Open questions and anything you could not resolve.
- PIDs of R jobs and exit status.
