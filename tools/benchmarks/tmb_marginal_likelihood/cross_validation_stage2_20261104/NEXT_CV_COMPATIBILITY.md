# Next CV compatibility work

Resume on inference-dev after the version 0.7.15 CV publication. Read RESUME.md,
full32-verified-state01.json and the current Git diff before changing source.
No owned numerical work remains. Do not replay closed producers or launch a new
panel merely because a saved script exists. Preserve other agents' work.

## 1. Compact training-selection snapshots

select_cv currently requires keep_fits = TRUE. Add an opt-in compact record of
native training log likelihood, complete K, used row IDs/outcomes/trials/weights,
formula/family/link, controls and diagnostics. Check these against retained native
fits before permitting selection without full models. Keep the current interface
and failure/tie policy. Prove identical choices and held-out scores while retaining
all candidates' warnings. This directly reduces memory without weakening guards.

## 2. Matched-set conditional-logistic / amt scoring

First qualify a narrow conditional-choice target using existing survival/amt
adapters. Keep every alternative from a choice stratum together in each split;
keep whole animals together when that is the intended transfer target. Refit
scales, coefficients and learned transforms inside analysis rows only. Score
normalized within-set choice probabilities with stable log-sum-exp arithmetic,
not unnormalized relative risk or independent Bernoulli probabilities. Match
native conditional likelihoods, weights, offsets, factor/interaction formulas,
ties and event units. Retain invalid strata and failures. Verify native K and
training-only selection on native oracles. Generic censored Cox prediction is a
separate target and must not inherit a conditional-choice scorer automatically.

## 3. Broaden ordinary-family and GAM targets

Qualify Gamma/NB GLM losses and their training-only dispersion/theta estimates
before extending Gaussian/Poisson/binomial scoring. Keep density versus mass,
constants, trials and weights explicit. For GAMs retain native formula metadata,
smoothing choices/effective degrees of freedom and multi-component family output
interpretation. Start with documented mean targets and native prediction oracles;
full predictive density requires its own qualified adapter.

## 4. Mixed and unmarked validation

Declare existing-group conditional versus new-group marginalized targets before
implementation. Setting random effects to zero is not integration. For unmarked
pcount, split sites/repeated visits according to the prediction question and
score the appropriate integrated observed-count likelihood with its used-visit
mask. Do not substitute latent abundance means or independent GLM visit scores.
Retain native nuisance/variance parameter counts and diagnostic gates.

## Acceptance and resource contracts

Current helpers construct folds; they do not qualify every model class for CV.
Keep the explicit split/response/trial/frequency-weight comparisons, all planned
assessment occurrences, retained failures/infinite losses and no post-hoc dropped
folds. Group-balanced score aggregation, inner-CV tuning and uncertainty/calibration
assessment require separately declared targets and validation designs.

Preserve fixed reference normalization, priors, interfaces, coefficient SEs, K,
all starts/cold/zero/limiting/curvature checks and private spatial/ML guards.
Use small declared native-oracle tests, generated help and focused feedback,
then a fresh normal vignette check and exact-archive installed suite for substantive
code changes. Source edits qualify old check coverage. Assertion totals are coverage
records, not a growth target. Reuse fixtures rather than adding duplicate fit panels.

Inspect actual R/source/resource identities before launch. Keep at most two owned
fit workers and three owned R interpreters including the driver, one thread each;
unknown inventory blocks launch, foreign jobs remain untouched and timing claims
remain ineligible during contention. Follow the disk policy: curated source only,
compact atomic records, no serialized full fits/chains/tapes or candidate arrays,
no broad cleanup. Retain failures and exact source/archive identities.

The previous operational controller expired at 07:00 Eastern November 4. Create
a fresh bounded controller for new work; do not alter or replay old producers.
The user is logging off now: this publication launches no further development.
