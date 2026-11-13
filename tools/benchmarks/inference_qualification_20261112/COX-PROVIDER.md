# Cox transformation metadata on Windows

The same four affected test files passed on R 4.6.2 when using survival
3.8-6. With survival 3.8-12, the existing 0.7.36 and 0.7.37 code reproduced
failures involving transformed retained Cox data. The newer provider
constructs an initial model frame with `na.pass` and applies its omission
function after evaluation. The matrix subset can discard basis classes and
attributes even when every row is complete. Retained prediction calls still
carry polynomial coefficients and spline knots in tested native fits.

The compatibility update preserves complete prediction calls, checks the
numeric basis after the provider discards its class, and uses the retained
fixed spline recipe for conditional-logistic designs. It checks complete
training rows and native design values before copying attributes from data
held during a package refit. Native coefficients, likelihoods, and
covariances keep their original checks.

For a native fit made outside the package, `na.action = stats::na.fail`
preserves the transformation attributes when the choice sets are complete.
A positive affine change to raw values can leave their standardized values
identical. A raw reference alone cannot recover the original center and
spread from those standardized values, so legacy predictions still need
the original recipe or a new fit with the original training data.

The provider-focused tests cover both R 4.6.1 with survival 3.8-6 and R 4.6.2
with survival 3.8-12. The R 4.6.1 focused suite passed 822 assertions. The
R 4.6.2 native Cox metadata test passed 157 assertions, and the broader
provider-focused tests passed on the same R version. Both comparisons used
the matched glmmTMB/TMB overlay and an explicit UTF-8 locale.
The R 4.6.2 follow-up candidate passed 43,162 assertions across all 240
installed test files and a separate 114 assertions in source-fixture tests.
The same archive rebuilt all vignettes in
`R CMD check --as-cran --no-tests --no-manual` and returned `Status: OK`.

A subsequent review identified a missing path for nested learned predictors
in the default refit backend. The targeted regression fits
`I(base::scale(z))` with `na.omit`, refits with a changed ecological
predictor, and independently reconstructs two uncentered predictions from
the original center and spread. The focused native tests passed 167
assertions on each provider. Unrecoverable nested recipes still use the
existing prediction error.

The nested-recovery candidate passed 43,164 assertions across all 240
installed test files, with the same 114 source-fixture assertions passing
separately. Its vignette-enabled package check also returned
`Status: OK`. The archive SHA-256 is
`039D66143CF9140C877C5E3BED6AD6D7F158906329240D214FB4161B6CD20A93`.
