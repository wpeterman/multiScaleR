# Selected native predictions, version0.7.9

Partial lme4 re.form predictions require fixed inputs plus grouping/slope inputs
for selected fitted effects, matching the native oracle across points, snapshots
and checked rasters. Native validation remains authoritative for fitted effects.
Contracts explicitly retain selected conditional effects and GAM terms/exclude/
iterms.type controls. Model penalties remain those of the full fitted model.
GAM excluded-term inputs remain required; its unchecked shortcut is unsupported.
No random-effect integration, kernel uncertainty or observation noise is invented.

Before the required metadata-only version bump: fresh package check Status: OK,
normal vignette/index rebuild, full exact-archive suite 6,924 passed / 0 failed,
12 unchanged warnings / 10 development-only skips. Focused:314 passed, eight
known development-load namespace warnings absent from the installed suite.
All 356 source hashes were verified. Only DESCRIPTION version and NEWS metadata
changed afterward; validation.json records both identities and qualification.
Priors, diagnostics and all starts/cold/Hessian/limiting gates are preserved.
Previous failures and unrelated local work remain retained separately.

Next compatibility checks should cover native GAM raster mean/component/SE
outputs with smoothing covariance and parallel partial-effect raster controls.
