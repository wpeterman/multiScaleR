# Fit-level diagnostics available from the archived study

The finite-grid study records four distinct sources of uncertainty.

1. The package can return no full bootstrap set if even one pair fails at
   any grid node. Always report the availability denominator alongside
   coverage conditional on return.
2. An available set can include truth while covering most of the search
   domain. Inspect its component count, bound contact, accepted width,
   and endpoint brackets before describing the scale as localized. A
   singleton contains only the accepted node. Its outside brackets are
   the adjacent rejected nodes and show the grid resolution.
3. Grid representation and bootstrap Monte Carlo error can affect the set
   independently. Nested coarse/fine comparisons hold the node decisions
   fixed. The 49/99 comparisons hold the dataset and node locations fixed
   and reuse the first 49 draws.
4. The scenario-level exact binomial intervals describe sampling across
   independent datasets. A within-dataset exact interval for a bootstrap
   tail describes simulation uncertainty at one node. Neither is a
   simultaneous confidence region for all grid nodes.

The root-refined chi-square set is used only to report overlap. It can
behave differently at weak effects and boundaries and supplies no
ground truth for bootstrap calibration.
