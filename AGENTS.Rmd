# AGENTS.md

## Purpose
This repository contains an R package. The goal of any automated coding work is to make targeted, reviewable improvements while preserving package stability, documentation integrity, and reproducibility.

## General rules
- Prefer small, minimal diffs.
- Do not make broad refactors unless explicitly requested.
- Preserve existing public function interfaces unless the task explicitly requires changing them.
- Do not guess about intended statistical or scientific behavior. Flag uncertainty clearly.
- Do not silently change numerical behavior, defaults, or return values.
- When touching package code, also inspect tests, documentation, and examples that may need updating.

## Package structure
Typical files and directories that may be relevant:
- `DESCRIPTION`
- `NAMESPACE`
- `R/`
- `man/`
- `tests/testthat/`
- `vignettes/`
- `src/` and `inst/` if present
- `.github/workflows/`

## Coding expectations
- Follow the existing code style in the repository.
- Keep functions focused and avoid unnecessary abstraction.
- Prefer readable base R unless the repository already uses tidyverse-style conventions in that part of the codebase.
- Avoid introducing new dependencies unless clearly justified.
- If a new dependency is necessary, explain why and update:
  - `DESCRIPTION`
  - namespace usage
  - any affected documentation/tests

## Documentation rules
- This package uses roxygen2 documentation unless the repository clearly indicates otherwise.
- If function behavior, arguments, defaults, return values, or examples change, update the roxygen comments in `R/`.
- Do not edit generated `.Rd` files directly unless explicitly instructed.
- Keep examples lightweight and CRAN-safe.

## Testing rules
- For any nontrivial code change, add or update tests in `tests/testthat/`.
- Prefer focused tests that verify behavior, edge cases, and regressions.
- Do not remove failing tests merely to make checks pass unless explicitly instructed and justified.
- If behavior changes intentionally, update tests to reflect the new contract and note that clearly.

## Validation expectations
After making changes, recommend the relevant validation commands. Common commands include:

```r
devtools::document()
devtools::test()
devtools::check()
```
