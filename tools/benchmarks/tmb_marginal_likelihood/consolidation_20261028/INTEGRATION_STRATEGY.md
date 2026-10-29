# One development line for scale inference

## Decision

Consolidate around the current master ancestry in a future `codex/inference-dev`
branch, with small reviewed commits for shared infrastructure and each method.
Keep the existing public entry points during the transition. This document is a
strategy; no branch switch, staging, commit or merge is authorized or performed.

The important gain now is a trustworthy choice of inference method for the same
model. More cache tuning has lower priority than independent posterior agreement,
clear uncertainty targets, and a common tested response likelihood.

## Branch reality

At inventory, master and `codex/bayesian-hierarchical-scales` both point to
87bf8779fe571921f2fd5d3e45c0f9e478530022. The performance and planar-kernel branches
also point there. The inference implementations are predominantly untracked
working-tree files, with supporting edits in DESCRIPTION, NAMESPACE, documentation
and existing interfaces. Git merge alone will not move this uncommitted work.

The old `dev` has 73 commits absent from this HEAD, while HEAD has 126 absent from
dev. `development` has 65 absent here and lacks 126 here. Do not merge either old
branch wholesale. Inventory any unique changes separately and port only changes
still needed, with tests. Do not reset or repurpose those branches.

Concurrent landscape-metric, native-prediction and other changes share this
checkout. Preserve their files and stage only explicit inference hunks when the
user eventually authorizes commits. The companion inventory records the exact
source identities used for this plan and flags shared files for manual selection.

## Share the model without conflating the inferential targets

For response data y, annual log scales a, population parameters theta and optional
response effects u, define one normalized conditional density p(y | a,u,theta).
The hierarchical likelihood integrates p(y | a,u,theta) p(a,u | theta) over latent
effects. A Bayesian posterior additionally uses explicitly declared proper
population priors p(theta), then integrates or samples theta as well. A Laplace
approximation describes an integration method; it does not determine whether the
target is likelihood or posterior.

| Layer | Shared contract | Method-specific behavior retained |
|---|---|---|
| Model specification | Response, trials, offsets, row identity, kernel design, fixed standardization, scale formula and grouping maps | Explicit capability checks and prior specification |
| Conditional response | Gaussian, Poisson, binomial, NB2 log density with all constants; NB2 size and Gaussian SD parameter meanings | Exact Gaussian response-effect marginalization uses a joint covariance adapter |
| Target | Immutable model identity and parameter map | Unpenalized population likelihood versus proper posterior; Jacobians belong to coordinates, not arbitrary penalties |
| Integration | Value, score, curvature and numerical accuracy metadata | Joint latent representation, Laplace, independent annual quadrature, exact Gaussian response integration |
| Solver | Start/seed/coordinate identity and retained failure records | Optimization, NUTS, importance correction, bootstrap |
| Uncertainty | Natural/log scale labels, conditioning and target population | Wald/profile confidence, repeated-latent prediction-error, posterior credible intervals |
| Diagnostics | No accepted result from an ineligible fit | Stationarity, curvature, limiting witnesses, quadrature error, importance ESS, NUTS diagnostics and functional precision |

Do not expose AIC for a posterior objective, call an empirical-Bayes plug-in
interval a full posterior interval, or use a positive equal-tail SD interval as
a test against SD=0. Keep data-only limiting witnesses and competing-optimum
rejections visible. Never select a method automatically because it returned a
number after another method failed.

## Current module map

| Component | Current source families | Consolidation action |
|---|---|---|
| Public Bayesian workflow | multiScale_bayes; bayes_scale_model/setup/methods/plot/simulate/prior_predictive/diagnostics | Preserve formula API and class; add explicit method metadata rather than rename now |
| Proper posterior computation | bayes_scale_likelihood/laplace/marginal/importance/beta_proposal/response_proposal | Reuse conditional response density; keep priors, Jacobians, population integration and correction diagnostics separate |
| Hierarchical likelihood | marginal_scale_model/likelihood/fit | Reuse model and response components, then remove preparation's dependence on Bayesian defaults in a later tested step |
| Accurate annual likelihood | factorized_scale_likelihood/fit/optimizer/score/limits/posterior/cache | Keep independent-year capability restrictions, exact zero-SD faces and strict cold audits; distinguish conditional annual summaries from full Bayes |
| Existing likelihood backends | profile_glm, joint_nb2, joint_mixed, conditional_logit, multiScale_optim, scale_inference | Preserve specialized profiled and conditional targets; connect through shared specifications only where identity tests establish equivalence |
| Kernel derivatives and preparation | kernel_derivatives, smooth_glm, kernel bins and supporting C++ | Later reuse after value/gradient/Hessian identity, pooled-standardization and tail-stability tests |

The factorized integrator currently requires one kernel and at most one scale
grouping term. Gaussian response intercepts must not couple scale groups.
Crossed annual effects, multiple kernels and additive groups can use the joint
Bayesian representation where supported, but do not inherit scalar-factorization
support. Temporal correlation, random slopes and custom families need new model
code and validation. Formula syntax does not create those capabilities by itself.

## Ordered integration batches and gates

1. Freeze source and input identities. Separate inference-only files from shared
   DESCRIPTION/NAMESPACE/interface/documentation hunks. Preserve research archives
   outside package builds and keep the original failures in every summary.
2. Share only the conditional response density first. Verify ordinary and AD
   values, scores, Hessians, normalization, offsets, trials, NB2 dispersion and
   posterior-versus-likelihood penalty differences. Keep stabilized kernel and
   Gaussian covariance adapters in place. Freeze a complete package check and
   the six scientific decisions before promoting the small batch.
3. Introduce a neutral internal model specification behind compatibility wrappers.
   Freeze the existing prepared-object contract, row maps and pooled centers.
   Keep saved-object compatibility and an explicit prior attachment step. Avoid
   copying Bayesian initialization rules into a population likelihood by accident.
4. Add explicit target/integration/uncertainty metadata and a capability table.
   Existing entry points continue working. Reject incompatible combinations
   before fitting. Use common diagnostics containers without hiding each
   method's original diagnostics or changing acceptance cutoffs.
5. Integrate research improvements independently. Exact response marginalization
   and transport require posterior density/Jacobian tests plus independent-seed
   agreement, not just better Rhat or zero divergences. Profile/limiting-model
   reliability needs cold score and witness checks. Each can ship separately
   when its gate passes; an unresolved method does not block unrelated fixes.
6. Expand predeclared same-input calibration across independent datasets, weak
   and strong signals, unequal years, zero/positive SD and valid prior
   sensitivities. Report cost, failures and Monte Carlo ambiguity, with
   dataset-level uncertainty. Only then choose stable public defaults or offer
   a unified convenience interface.

Each batch should be independently reviewable and reversible. When commits are
authorized, bump the package version once for the related integration batch,
update NEWS, and attach exact validation evidence. Do not mechanically commit
the whole dirty checkout or copy whole shared files over concurrent work.

The five-file response batch depends on the larger uncommitted inference
foundation. It is not a standalone cherry-pick onto HEAD. Before constructing
the first integration commit, include the dependency-complete model preparation,
kernel engine and native support it actually calls. RcppExports and native
registration must match the included C++ files. Check that exact integration
tree again, since a check of the full dirty source does not validate an arbitrary
selected subset. Separate managed worktrees can keep subsequent concurrent
feature edits from changing an integration candidate during its gate.

Expose chain-based probability MCSE for a user's specified scale thresholds in
a later diagnostics batch. Retain current passed flags and unavailable reports
for compact/approximate fits. Constant indicators need a nonestimable result,
and probability precision must be judged against the scientific decision.
The long-reference results show why raw Rhat/ESS alone are insufficient.

## Concurrent prediction overlay discovered during this study

Independent work under tools/benchmarks/bayesian_scales/prediction_20261027
adds predict_scales() and new-observation Bayesian prediction/simulation, with
explicit new-group policies and frozen training designs. Its own report and
validation must determine readiness. This study did not implement that feature.
Its latest source is included in the final package check, while its prospective
predictive calibration has a separate owner and evidence. Include the additional
source, generated help, tests and prediction
vignette in the eventual integration inventory; concurrent-source-delta.json
records the first observed identities against this study's original snapshot.

Preserve those changes when consolidating. In particular, do not overwrite the
new prediction metadata or current public methods with older frozen files merely
because those files passed earlier checks. The shared response layer's live
promotion requires a fresh source gate covering the relevant current code.
That full check now passes, with strict source equivalence supporting reuse of
the completed scientific evidence. Existing numerical evidence can be reused only after explicit
source/input equivalence; source drift alone is not a reason to rerun identical
numerical fits. Full package checks need a free interpreter window and must not
disrupt the other owner's predeclared example/predictive workers.

The five response files are now integrated locally and pass 1,316 live assertions.
The final audit verifies their checked source and all protected code/help/tests.
A descriptive concurrent NEWS addition is the sole inspected live-source
exception. Operational isolation remains a separate unfinished gate: the original
guard missed three existing rsession-utf8.exe hosts. Its corrected name and
unavailable-inventory tests pass, and it now blocks a new full check while those
sessions remain. Do not stop foreign RStudio sessions to satisfy this gate.
Use an actually free window for the eventual dependency-complete integration
check; do not relabel the earlier check as isolated or repeat its completed
numerical fits.

## Where continued work has scientific value

The exact transport substantially improved eligibility in the ten stronger
holdouts (10/10 versus 4/10 for the original runs), but selected eligibility is
not general calibration. All four long weak-case references now pass their
original and functional precision gates, with no flagged annual/SD endpoint or
median disagreements across four independent comparisons. This clears the
selected agreement gate. Public coordinate reconstruction and broader independent
calibration remain separate requirements.

Pooling helps when years share information and individual annual scales are
weakly measured. Both hierarchical likelihood and Bayes can pool. Full Bayes
also propagates population-SD uncertainty and allows explicit prior information;
those benefits can matter with few years or weak signals. Priors can also hurt
point recovery. The four weak unequal-year illustrations favored unpooled point
estimates on average, so they do not establish a universal Bayesian advantage.

Stop polishing low-level speed once it no longer changes feasible models,
reliable acceptance, or useful uncertainty. The next substantive milestones are
posterior agreement, boundary-aware SD uncertainty, independent calibration and
one understandable workflow with honest method limits.
