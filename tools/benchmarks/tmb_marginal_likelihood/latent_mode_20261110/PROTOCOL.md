# November10 release and latent-mode qualification

User authorized commit/push to inference-dev, followed by GDR latent-mode
qualification, preserving priors, diagnostics, parameter counts, guards and
uncertainty distinctions. Logoff17:30 Eastern. Stop long launches17:00, all
numerical launches17:10, owned work17:15, checkpoint17:20. Use control01.json,
one thread, maximum two fit workers/three owned R, known foreign/RAM exceptions,
complete live inventory and disk reserves. Leave foreign jobs untouched.

First publish the already-validated development batch as0.7.32. The patch-level
version bump follows the repository rule after code validation. Shipped R code,
namespace, tests, generated help and vignette sources must match the full29
archive (42,575 passes, ten unchanged skips, zero findings). Build and normally
check the new version with vignettes, proving archive code/test identity against
full29. Its full test coverage is inherited by exact code identity; do not claim
that the entire installed suite was rerun on the metadata-only version bump.
Exclude the owner-edited unified DEVELOPMENT.md and unrelated artifacts from
staging. Push inference-dev only; no merge or tag.

After push, declare the fixed-point latent-mode audit in a separate protocol
before launch. The latest memo is preserved in memo-reviewed.md with hash.
Prioritize retained competing inner branches and direct current joint values/
gradients/curvature. Constructor seed zero under a cached inner policy differs
from forcing zero at every evaluation. No automatic rescue, passing-seed search,
prior relaxation, covariance repair, sampler/profile/recovery panel or public
range/CV promotion. Preserve every failed/interrupted result and checkpoint only
compact scalar diagnostics and identities, never full models/tapes/node arrays.
