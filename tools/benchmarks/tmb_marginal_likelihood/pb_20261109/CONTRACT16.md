# Private prepared-binomial prerequisites in 0.7.29

This internal context consumes ordered raw numeric features from a retained
ordinary logit-binomial GLM. Literal physical spline knots/boundaries, factor
coding, products and fixed offsets stay at their training definitions. The
context retains numeric recipes and target identities, not the original fit,
caller formula environment or mutable coefficient cache.

Binary numeric/logical/two-level-factor responses and two-column integer
success/failure responses use positive integer frequency weights. The target
is the sum of frequency-weighted, normalized binomial log probabilities.
Proportion/trial-weight conventions and zero/noninteger weights are not
silently reinterpreted by this private adapter. Reference reconstruction and
exact ordered row/column/target identity are mandatory. A declared source hash
does not authenticate external feature geometry.

Coefficient existence is separate from IRLS convergence. Exact rational checks
use the original represented double design. A full-rank paired outcome basis
can establish overlap without an LP. Otherwise lpSolve proposes a balancing or
recession witness and optional gmp verifies it exactly, including rank and
strict positivity. Structural zero pivots use explicit exact row pivoting.
A non-strict sum maximizer does not establish quasi-complete separation; the
subtype needs its own proof. Missing arithmetic, unstable proposal bases and
unverified witnesses remain unknown. Every solver failure is retained.

Accepted coefficient points also require finite normalized likelihood, full
native numerical rank, convergence, no fit warnings and the unchanged scaled
coefficient-score gate. They grant no outer-optimizer, model-selection or
uncertainty admission. Nominal counts retain all coefficient columns plus all
declared scale coordinates, including the five-plus-two and seven-plus-three
examples. No coefficient or scale SE is manufactured by this value-only packet.

The raw analytic, hard-radius, authenticated-provider, AIC/CV and uncertainty
guards remain in force. Existing public models, priors, normalized targets,
fixed reference transforms, diagnostics and parameter penalties are unchanged.
Next work is authenticated providers and complete outer-search qualification,
then downstream prediction/selection/CV integration and equal-work benchmarks.
Spatial response fields are separately prioritized in SPATIAL_PLAN.md.

Qualification: 252 installed focused expectations and fresh-worker equality;
32 prespecified geometry cases agreeing with an independent exact rcdd LP;
fresh normal vignette-enabled check and complete exact-archive installed suite.
See VALIDATION16.md. These are software and fixed-point contracts, not empirical
recovery or optimizer-efficiency claims.
