# Validated 0.7.32 integration

This patch integrates training-only Gaussian spatial-field CV with estimated
range and ecological scales, joint held-out scoring and guarded model selection.
It adds private normalized GDR range/CDF prerequisites, optional retained inner
controls, typed optimizer checkpoints before covariance and current local-mode
diagnostics. Public estimated GDR range/CV and unqualified inference remain
guarded. Priors, parameter penalties and native uncertainty scopes are preserved.

The statistical source passed a normal vignette-enabled check and42,575
exact-archive assertions, with zero failures/errors/warnings and ten unchanged
skips. Parent archive SHA256:
0980a080a27787d5fa1618604a1c7509126ba1efa2bf27fb4c4c560801cf6419.
The feedback planner passed12 Python tests and retains every test in full mode.
The independent small joint-tape audit matched native values/gradients exactly
and random Hessians within3.41e-13. Original failures remain retained.

Version0.7.32 changes release metadata after statistical code qualification,
following the repository's patch-version rule. Its normal vignette-enabled
check isOK. Every shipped R/namespace/help/test/native/vignette source member
matches the fully tested parent archive. The entire installed suite was not
rerun for the version-only change. Fresh release archive SHA256:
69195a82d4c470f4d08beee4e5f7a4e52527402135748b5ac644ecbde2aa78d9.

The first release driver stopped after a CRLF-only version-field regex failure;
the completed build was preserved and its first normal check continued under
release05. No statistical code, test or numerical gate changed.

The user's explicit commit/push instruction covers this validated batch.
Unrelated unified DEVELOPMENT.md edits are excluded. NEXT_PRIORITIES.md orders
the latest GDR memo requirements. Mode qualification after publication remains
separate from new fitted-model or density admission; no automatic rescue is
claimed. Operational data, model/tape arrays and large validation artifacts
are kept outside this commit.
