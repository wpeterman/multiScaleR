# Private spatial source preflight

The declared optional report is closed and validated at02:07 Eastern on October29. The final317-file candidate passes its complete gate below. Earlier progress notes remain as history.
It reads supplied CRS, grid geometry and IDs without reading raster values or
re-extracting landscapes. Coordinate comparison uses explicit keys and a recorded
absolute tolerance (exact zero by default), without transformations or rescaling.
Nominal buffer/extent geometry remains distinct from actual valid-cell/kernel
coverage. Historical source and extraction-unit lineage cannot be certified from
compact inputs. No inference eligibility, fitting, priors or uncertainty changes.

unit01 has34passing assertions and one error in the missing-raster-CRS rejection
path. sf raised its own empty-CRS parsing error rather than the helper's intended
message. That failure and original source/tests/logs are preserved. The helper
now catches that parse failure only to reject the source explicitly. New tests
cover the same declared fixture and additional malformed geometry/tolerance/ID
cases. No source or numerical acceptance gate in an existing method changed.

The prior 5,261-assertion clean check belongs to response_sd_profile_20261029's
frozen314-file candidate. It does not cover this new R/test/help source or NEWS
edit. Fresh generated help, compatibility and complete candidate validation are
required. No new sampler, response fit or scientific dataset is declared.

At01:47 Eastern, unit02 passes64assertions with zero failures, errors, warnings
or skips. Missing-CRS rejection is explicit, and all declared source/row/coordinate,
geometry/tolerance and nominal extent cases pass. The raster-hole equality test
confirms no valid-cell/kernel-coverage inference follows from grid metadata.
Documentation, compatibility and fresh full validation remain pending.

At01:50 Eastern, generated internal Rd passes checkRd and NAMESPACE is unchanged.
compatibility01 passes 263 assertions across six existing spatial,
preparation, grouped-raster and prediction-contract modules, with zero failures,
errors, warnings or skips. unit02 remains64passes. A fresh317-file candidate is
being frozen; only NEWS changes among prior files, and helper/test/generated Rd
are new. The full candidate check and exact installed tests are still pending.


## Final validation

The fresh317-file source candidate passes the standard as-cran check, including
examples, donttest examples and rebuilt vignettes. Its exact unchanged archive
tests ran separately against the check-installed package with NOT_CRAN=true:
5,325passing assertions, zero failures,12existing warnings and10development skips.
No new test is skipped. Focused64 and compatibility263 are included in that total.
The full candidate gate has45resource samples, maximum3actual R interpreters,
minimum13.88GiB RAM and46.22GiB free disk, with no unknown or foreign processes.
No R process remains at closeout.

Archive SHA256:
118A69E2A54465CD56EE833A85F80C0B3A379BA0D7AA25084190CC3B40A3FC8A.
All317source identities and extracted archive-test files are exact. Branch, HEAD,
version and empty index remain unchanged. All five closed predecessor manifests
and historic pinned evidence pass. NEWS is the only modified old source; the
private helper, test and generated help are new. Every existing fitter, prior,
acceptance/diagnostic gate, prediction value/attribute and public interface is
unchanged. The original missing-CRS failure/source remain in this evidence set.

The report verifies supplied current metadata and coordinate alignment. It does
not certify historical source/unit lineage from compact inputs, evaluate raster
holes or actual weighted support, reject a fitted model on nominal buffer
clipping, fit either target or establish uncertainty calibration. The next guide
is separately declared; its latest plan uses the existing bundled example and
stored fit as read-only consumers. No closed source or producer is rewritten.
