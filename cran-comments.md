## Test environments

Windows 11 x64, R 4.6.1: the final 0.7.3 source archive passed
`R CMD check --as-cran --no-manual` with 0 errors, 0 warnings, and 1 note.
All examples, installed tests, and six vignette rebuilds passed.
The only note was "Days since last update: 4"; the reason for this
corrective release is explained below. Full development tests also passed.
Other platforms and the PDF manual have not been rechecked for this update.

## Changes since the previous CRAN release

Version 0.7.3 corrects parameter counts used in AIC, AICc, BIC, and sigma
profiles for unmarked models. The previous fallback could count predictor
names rather than all fitted parameters. For a reproduced random-intercept
count model, the base count was incorrectly 3 rather than 8, and the total
including one estimated scale was 4 rather than 9. This changes information
criteria and can affect model rankings, but does not change fitted
likelihoods, coefficients, or scale estimates.

Native unmarked likelihood and parameter metadata are now validated, with
focused tests for multiple model families, mixtures, fixed parameters, and
random effects. Undefined AICc values remain NA in profiles; model-selection
tables provide an explicit error for insufficient sample sizes.

This corrective release is needed shortly after 0.7.2 because the previous
parameter counts can give misleading model-selection results.

Examples and vignettes use bundled data and cached fitted examples without
network access. Private data used to confirm the repair are not included.
