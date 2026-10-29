# Private spatial source and coordinate preflight

Declared before implementation or new numerical work. Start only after the
response_sd_profile_20261029 source/evidence batch is closed. Remain on
codex/inference-dev without staging, commit, bump, merge or push.

## Why this is next

The private spatial basis accepts two numeric columns described as projected
coordinates, but a plain data frame cannot establish their CRS or linear unit.
The compact kernel inputs retain point IDs, distance scaling and binned values;
they do not retain original point geometry, raster CRS or verified historical
source lineage. Numeric coordinate magnitudes cannot supply that information.
Do not infer longitude/latitude from ranges or treat a CRS label as verification
of the original raster values. This affects both likelihood and posterior spatial
preparation and can be checked without fitting or additional scientific datasets.

## Bounded implementation

Add one optional private source-report helper, with explicit original point and
raster objects, compact inputs, observation data and coordinate/row-ID names.
Keep all existing public and private fitting/preparation/prediction interfaces,
numerical outputs, priors, gates and frozen source studies unchanged. Do not
automatically call the helper from an existing fitter or attach new attributes
to existing prediction outputs.

Validate defined projected POINT geometry and a defined matching raster CRS,
report the authoritative CRS linear-unit label, and align source point IDs to
the binned IDs and observation rows by key. Repeated observation IDs are allowed
when they refer to the same source point. Reject missing/duplicate source keys,
unmatched rows, nonfinite coordinate columns and explicit coordinate disagreement.
Compare the supplied columns with original point coordinates without transforming,
rescaling, snapping or remapping. Default coordinate tolerance is exact zero;
any caller-supplied nonnegative absolute tolerance is recorded in map units.
Support sf and terra point objects through their ordinary geometry adapters.

Return compact source/grid/unit/alignment summaries and per-source-point nominal
buffer/extent diagnostics at inputs$max_D. Distinguish an extent-interior buffer
from actual valid-cell or kernel-weight coverage: raster holes, masking and
covariate-specific support remain unknown unless separately evaluated. A clipped
nominal buffer is descriptive evidence, not automatic fit rejection. Do not read
full raster values, copy binned/candidate arrays, re-extract landscapes or change
their weights. The report must explicitly say that historical raster-value
lineage is unavailable from compact inputs and is not certified by this check.
Compact inputs also do not prove the original extraction's distance unit. Any
nominal buffer diagnostic must record its conditional assumption that max_D is
expressed in the supplied projected source's linear unit; it cannot certify
that historical unit agreement. Do not relabel the internal scaling factor as
a conversion between meters and feet. An explicit future provenance capture
would need a separate small preparation batch, preserving legacy objects.

Keep field covariance range, ecological kernel scale and inputs$unit_conv
(internal distance scaling) distinct. Do not infer scientific support, inference
eligibility, calibrated uncertainty, accepted nonstationarity or a fitted model
from the report. Population priors and target-specific diagnostics are unaffected.
Add generated internal help and a brief NEWS entry; no new dependency or export.

## Focused fixture and required assertions

Use one deterministic small projected raster with two named source points and
matching compact Gaussian kernel inputs. No response fitting, sampler, seed
search, scientific replication or timing claim is declared. Check valid sf and
terra sources, explicit keys/reordering/repeated observations, matching coordinate
columns, retained inputs/RNG, missing/geographic/different CRS, non-point geometry,
invalid IDs/columns/tolerance, clear coordinate-unit mismatch and extent-edge
reporting. Include a raster-hole example proving that extent-interior status
does not assert valid-cell coverage. Assertions concern software contracts only.

Before numerical launch, inspect known source/process/resource identities. One
owned worker/thread, max four actual R interpreters globally, RAM >12 GiB and
disk >5 GiB; unknown inventory blocks launch. Bound focused phases to600seconds.
Preserve all failures and compact atomic records. Then generate internal Rd,
run focused spatial/preparation compatibility checks and a fresh exact full
candidate gate in an R-free window. Do not reuse a closed check for changed
source or alter closed baselines/manifests. No long launches after04:00 Eastern,
first wake at/after04:15 stops launches, safe checkpoint by04:35 and delete the
one-night heartbeat by04:45. After cutoff checkpoint/delete only.
