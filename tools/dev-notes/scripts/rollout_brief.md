# Vignette rollout brief (multiScaleR, branch docs-restructure)

You are rewriting one or more vignettes of the R package `multiScaleR` for its
maintainer. He reviewed and approved two pilot rewrites. Your job is to apply
the same pattern to the vignette(s) named in your task. He will review a
rendered HTML of each before anything is committed.

## Where things are

- Worktree (work ONLY here, never `cd` to the parent repository, never switch
  branches): `C:\Users\peterman.73\OneDrive - The Ohio State University\R\Packages\multiScaleR\.claude\worktrees\docs-restructure`
- Scratchpad: `C:\Users\peterman.73\AppData\Local\Temp\claude\9e3b7d3898bd820b\459cb932-8333-43df-b0da-7d9a2422aa59\scratchpad` (called SCR below)
- The approved pilots, which are your pattern. Read both fully before writing:
  `vignettes/quickstart.Rmd` and `vignettes/inference-contracts.Rmd`.
- Originals of every vignette: `git show bfdc133e3:vignettes/<file>`.
- `tools/dev-notes/vignette-relocation-map.md` lists blocks removed from
  `inference-contracts.Rmd` and where they belong. If a row names your
  vignette as destination, bring that content in (rewritten in the pattern),
  reading the original text with `git show bfdc133e3:vignettes/inference-contracts.Rmd`.
  Do NOT edit the map; report which rows you completed.
- `tools/dev-notes/validation-status.md` is a maintainer note. Development
  status narrative (validation gates, private references, "not yet
  calibrated" histories, defect reports about other packages) leaves the
  vignette. APPEND such passages verbatim to that file under a heading
  `## From <vignette file>` with the original line numbers. Keep in the
  vignette only what a user needs: a short "Current restrictions" table or a
  one-sentence statement that a feature is experimental.

## House style (strict; the maintainer checks this)

Read these skill files first and follow them:
`C:\Users\peterman.73\claude_skills\peterman-research-suite\skills\academic-voice-style\SKILL.md`,
`...\skills\r-doc-style\SKILL.md`, `...\skills\r-vignette-writing\SKILL.md`,
`...\skills\humanizer\SKILL.md` (mandatory final pass).

Non-negotiable mechanics:
- No em-dashes anywhere. No semicolons in prose. Colons sparingly (the bold
  "**How to read this output:**" labels are the accepted exception).
- American English, active voice, second person ("you"), concise.
- Never say "cell" for a simulation design combination; say "scenario".
  (Raster cells are still cells.)
- Headings: `#` top level, sentence-case task phrases ("Prepare the data",
  "Read the summary"). Not noun stacks.
- Spell out Greek letters in prose and tables (sigma, beta, tau). Equations
  keep notation.
- Define each technical term inline on first use, in plain words.
- After every printed output or figure that a reader must interpret, add a
  short "**How to read this ...:**" block that says what the numbers mean for
  the study, in the units of the data. Translate estimates (for example, an
  effective distance in meters, a multiplicative change in expected count).
- Blockquote callouts for consequences: `> **What this means for your study.** ...`
  Use them sparingly (one to three per vignette).
- Every code chunk is annotated with short comments. Visible code must be the
  real call a user would type. No helper functions such as `pkg_extdata()` in
  visible code; use `system.file("extdata", ..., package = "multiScaleR")`.
- `terra` first for spatial work. Use `sf` only where a function requires it,
  and say so.
- Do NOT invent or add literature references. Keep existing citations exactly
  as they are (they were verified). Do not fabricate claims about methods. If
  you are unsure what a function does, run it or read `R/` and `man/`.
- Do not state anything about package behavior you have not verified by
  running code or reading the source. Remove rather than guess.
- Marginal-effect plots no longer print a caption by default
  (`plot_marginal_effects(fit, add_label = TRUE)` adds it). Kernel plots
  print three lines under each panel by default; `plot(fit, add_label = FALSE)`
  removes them. Do not describe a note under marginal-effect panels.

## Structure every vignette should have

1. `# What this vignette covers` (two short paragraphs: the question it
   answers, what you will do, what you need to have read).
2. Task-ordered sections with worked, runnable examples and interpretation.
3. `# Quick reference`: one code block that matches the body exactly (same
   object names, arguments, seeds), then `## Summary of key functions` (a
   table), then `## See also` (bulleted links to other vignettes written as
   `vignette("file-name", package = "multiScaleR")` using CURRENT file names;
   files are renamed in a later phase).
4. YAML as in the pilots (`rmarkdown::html_vignette`, `number_sections: true`,
   `toc: true`), title equal to `\VignetteIndexEntry`. Keep the current file
   name. You may improve the title to a plain task phrase.

## Examples and data

- Use the running example wherever the feature allows: `data("landscape_counts")`,
  `data("surv_pts")`, `terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))`.
  Truth: sigma 250 m for `land1`, 500 m for `land2`; effects land1 = -0.5,
  land2 = 0.7, site = 0.3, intercept 0.25; `land3` has no effect.
  `surv_pts` carries no CRS; functions that need a projected CRS take `crs =`.
- Fits that run in a few seconds run live (`backend = "profile_glm"` with
  `kernel_prep(bin = TRUE)` fits the running example in about 3 s; the default
  `refit` backend with binned inputs is also fast). Target under 60 s total
  build time per vignette.
- Slow fits: show the real call with `eval = FALSE`, load the stored result
  in a hidden chunk. Reuse existing objects in `inst/extdata/` (see
  `data-raw/generate_vignette_cache.R` for what `vignette_cache.rds` holds).
  Do NOT edit `inst/extdata/vignette_cache.rds` or
  `data-raw/generate_vignette_cache.R` (other agents share them). If you truly
  need a new stored object, write `data-raw/cache_<vignette>.R` that creates
  `inst/extdata/<vignette>-cache.rds` (keep it small, xz compressed), run it,
  and report its size.
- No global `eval = FALSE`. No chunks that silently print nothing. Every
  shown result must come from code that ran or from a stored object.
- Load stored objects with the pilot's pattern (local path first, then
  `system.file`), in a hidden setup chunk.

## Hard limits

- Edit only your assigned vignette file(s), `tools/dev-notes/validation-status.md`
  (append only), and any new `data-raw/cache_*.R` / `inst/extdata/*-cache.rds`
  you create. Do not edit `R/`, `man/`, `tests/`, `DESCRIPTION`, `NAMESPACE`,
  `NEWS`, `_pkgdown.yml`, the README, or other vignettes.
- Do not commit, stash, switch branches, or push. No `git` write commands.
- Do not install packages. Do not modify the R library at `SCR/rlib`.
- Job discipline: launch each R job so you know its PID
  (`cmd > log 2>&1 & echo "PID $!"; wait $!`), and only ever wait on or stop
  that exact PID. Never use `pkill`, `killall`, `taskkill /IM`, or
  `Stop-Process -Name`. Never touch a process you did not start.
- No user name, email, or project details in any request to an outside service.

## How to render and check (Git Bash)

```bash
cd "<worktree>"
SCR="C:/Users/peterman.73/AppData/Local/Temp/claude/9e3b7d3898bd820b/459cb932-8333-43df-b0da-7d9a2422aa59/scratchpad"
RS="/c/Program Files/R/R-4.6.2/bin/Rscript.exe"
# render against the private build of this branch (absolute output dir required)
"$RS" --vanilla "$SCR/render_one.R" <file>.Rmd "$SCR/render"
# plain text + figures for inspection (view PNGs with the Read tool)
"$RS" --vanilla "$SCR/html_inspect.R" "$SCR/render/<file>.html" "$SCR/render/<file>"
# style audit on prose only
"$RS" --vanilla "$SCR/prose_only.R" vignettes/<file>.Rmd "$SCR/render/<file>_prose.md"
"/c/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe" -B C:/Users/peterman.73/claude_skills/peterman-research-suite/skills/humanizer/scripts/humanizer_audit.py "$SCR/render/<file>_prose.md"
```

For exploratory R scripts, write them to `SCR/<your-vignette>_explore.R` and
start with `.libPaths(c(file.path(scr, "rlib"), .libPaths()))` so the branch
build (0.7.36 plus the caption change) is used, not the old user-library
version. Windows shell quoting breaks `Rscript -e` with pipes, so use script
files.

Read the rendered text and look at every figure. Check each interpretive
sentence against the numbers actually printed. Fix anything that does not
match. When the render is clean, copy the HTML to
`<worktree>/.claude/pilot-review/<file>.html`.

## What to report back (concise)

- Line counts and prose word counts, original versus new.
- Section outline of the new vignette.
- What moved out and where (dev note, another vignette's owner, a help page),
  with original line ranges, so the relocation ledger can be updated.
  Content that belongs on a help page: do not write it into `R/`; list it with
  original line ranges and the target help topic.
- Relocation-map rows completed.
- Every claim you could not verify and removed, and anything in the original
  that looks wrong (code that errors, text contradicting behavior).
- Render time, any warnings, humanizer audit result.
- PIDs of R jobs you launched and their exit status.
