# Refit reporting fixes, version 0.7.7

The refit profile cache uses bounded data-specific keys and validates its payload.
Invalid active Hessian curvature yields missing SEs with retained diagnostics and
raw Hessians. Quiet fitting suppresses printed diagnostic blocks. Selection
labels separate kernel from response family. Native likelihood degrees of freedom,
including dispersion and GAM effective degrees of freedom, determine penalties;
active spatial parameters count even when their SE is missing.

Validation before the required metadata-only version bump: fresh vignette-enabled
package check Status: OK, full exact-archive suite 6,846 passed / 0 failed,
12 unchanged warnings / 10 development-only skips. Focused tests: 307 passed,
no failures/warnings/skips. Generated documentation is refreshed. Only the
DESCRIPTION version and NEWS section changed afterward; code, help, tests and
examples retain their exact validated hashes. validation.json records both
identities. This establishes software behavior without new calibration claims.

All intermediate failures and unrelated local work are retained separately.
Priors, numerical fit gates and uncertainty distinctions remain intact.
Next priorities: native zero-random-effect row requirements, followed by
provider-specific GAM snapshot controls with direct/prepared/raster parity.
