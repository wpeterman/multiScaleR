# Landscape-metric expansion roadmap (`landscape_var`)

Scoping memo for adding categorical landscape metrics to `multiScaleR`, parallel to
the continuous surface-metric work in `surface_var()`. Validation scripts that back
this memo live in `bench_tmp/` (`.Rbuildignore`d, so excluded from the package build):
`L1_reference.R` (formula validation), `L2_edge.R` (edge behaviour), `L3_fft.R` (FFT
projection proof + benchmark), `L4_clumpy_closedform.R` (CLUMPY closed form),
`L5`/`L6` (the implemented direct and FFT paths vs `landscapemetrics`), with shared
reference implementations in `Lref.R`.

**Status (v0.6.36): the priority-1 through priority-4 additions below are IMPLEMENTED**
and shipped: `iji`; the information-theory family `ent`/`condent`/`joinent`/`mutinf`/
`relmutinf` (with the `base` argument; `base = 2` for bits); and the class-level
`clumpy`/`pland`/`ca` via a new `class =` argument to `landscape_var()`. All reproduce
`landscapemetrics` to machine precision through both the fixed-buffer C++ path and the
FFT raster projection. Remaining out-of-scope items (edge contrast, nLSI) are noted
below.

Reference: Cushman, McGarigal & Neel (2008) *Ecological Indicators* 8:691–703,
"Parsimony in landscape metrics"; Nowosad & Stepinski (2019) *Landscape Ecology*
34:2091–2101 (information theory). Validation target: `landscapemetrics` 2.2.1.

---

## 1. The gate: FFT-projectability

`kernel_scale.raster` projects a metric to every cell via FFT convolution only if the
metric is a function, **per moving window**, of one of:

1. **composition** — per-class cell counts (`.landscape_composition_raster_fft`)
2. **edge** — edge-segment counts (`.landscape_edge_counts_raster_fft`)
3. **configuration** — the class-pair adjacency / co-occurrence matrix
   (`.landscape_adjacency_raster_fft` → `.landscape_pair_count_raster_fft`)

Anything that needs connected-component **patch labeling**, a **distance-transform**
core, or **nearest-neighbor** search is NOT projectable by this machinery.

The adjacency pipeline already builds the full per-window class-pair co-occurrence
matrix (it powers `contag`/`ai`/`pladj`), so any metric that is a closed-form
function of that matrix is essentially **free to add** — only O(cells) arithmetic on
top of the convolutions already being run.

**Current set (16):** composition `shdi shei sidi siei msidi msiei pr prd rpr ta`;
edge `ed te lsi`; adjacency `ai pladj contag`.

---

## 2. Recommended additions (priority order)

| # | Metric(s) | Level | Reuses | Cushman / rationale |
|---|-----------|-------|--------|---------------------|
| 1 | **IJI** | landscape | adjacency | Interspersion & juxtaposition; a near-universal landscape configuration component. Pure function of the between-class co-occurrence matrix. |
| 2 | **ent, condent, joinent, mutinf, relmutinf** | landscape | adjacency | Nowosad & Stepinski information theory; cleanly separates composition (marginal entropy) from configuration (mutual information). All from the normalized co-occurrence matrix. |
| 3 | **CLUMPY** | class | adjacency + composition | Composition-corrected aggregation, bounded [-1, 1]. Needs a new `class=` argument (class-specific). |
| 4 | **PLAND / CA** | class | composition | Basic composition; partly redundant with `kernel_var()` on a 0/1 layer but convenient and standard. Needs `class=`. |
| 5 | nLSI | landscape | edge | Minor; normalized LSI. |

**Edge contrast (CWED / TECI / ECON):** Cushman's single most universal *class* metric
and FFT-able (weighted edge sums), but requires a user-supplied class-pair contrast-
weight matrix — new infrastructure, and not present in `landscapemetrics` (FRAGSTATS-
only). High ecological value, higher implementation cost. Deferred.

**NOT projectable (patch-based; Cushman-important but out of scope):** NP, PD, LPI,
AREA_*, MESH, DIVISION, SPLIT, COHESION; SHAPE, FRAC, PARA, CONTIG, CIRCLE, PAFRAC,
GYRATE; ENN, PROX, SIMI; patch-level core CORE_MN/CV, DCORE, NDCA, DCAD, CAI_MN/CV.
(TCA/CPLAND/CAI aggregate-core could be approximated only via depth-1 erosion; the
patch-level core distribution is not projectable.)

---

## 3. Validated formulas (exact conventions)

All reference implementations reproduce `landscapemetrics` 2.2.1 **exactly**
(`L1_reference.R`: max abs diff = 0 over 5 random 30×30 four-class rasters for
`iji ent condent joinent mutinf relmutinf clumpy pland`). The conventions that matter:

- **Adjacency:** double-counted rook (4-neighbour) adjacency (FRAGSTATS convention).
  `A[i,k]` = number of physical i–k edges; symmetric.
- **IJI:** `m` = number of classes *present in the window*. `NA` when `m < 3` or when
  there are no between-class edges (`E = 0`).
  `IJI = -Σ(p·ln p) / ln(0.5·m·(m-1)) · 100`, `p = e_ik / E` over unordered pairs.
- **Information theory:** co-occurrence `P = A / sum(A)`, log base 2.
  `ent = H(rowSums P)`, `joinent = H(vec P)`, `condent = joinent − ent`,
  `mutinf = ent − condent`, and crucially
  **`relmutinf = 1` when `mutinf == 0`** (else `mutinf / ent`) — the `landscapemetrics`
  `ifelse(aggr==0, 1, aggr/comp)` convention.
- **CLUMPY:** pad the whole matrix with a `-999` background border, double-count
  adjacency, `min_e` = minimum **perimeter** (`4n` / `4n+2` / `4n+4`, not max-like-
  adjacency), `g_i = like_i / (other_i − min_e_i)`, piecewise CLUMPY. `NA` when the
  class fills the window (`p_i == 1`) or is absent.
- **PLAND:** `100 · n_class / Σ n` (NA cells excluded from the denominator).

These exact conventions are pinned in `bench_tmp/Lref.R` and should be ported directly.

---

## 4. Edge behaviour (`L2_edge.R`)

Scenarios: single-class (m=1), two-class stripes, three-class random, interior-NA
holes, rare single-cell class. Reference matches `landscapemetrics` in **every case
except one**:

| Scenario | Behaviour (all match `landscapemetrics` unless noted) |
|----------|-------------------------------------------------------|
| single class (m=1) | IJI = NA; ent = mutinf = 0; **relmutinf = 1**; clumpy = NA; pland = 100 |
| two-class stripes | IJI = NA (m=2 < 3); ent = 1; mutinf = 0; relmutinf = 1; clumpy, pland defined |
| three-class random | all defined and exact |
| rare single-cell class | all defined and exact (clumpy = 1 for the isolated cell) |
| **interior-NA holes** | all exact **except CLUMPY** (ref 0.1090 vs lm 0.0952) |

**The one nuance — CLUMPY with interior NA:** the formula is exact on complete rasters;
the disagreement is purely how interior `NA` cells enter the co-occurrence (treated as
background vs excluded). When implemented over the FFT pipeline, CLUMPY's NA handling
will follow the package's existing `na.rm` / `.landscape_clean_count_matrix`
convention and be re-validated against `landscapemetrics` at that point. Invariant
landscapes (m<3, single class) already return NA/identity consistently.

---

## 5. FFT projection proof + benchmark (`L3_fft.R`)

**(A) Projection is exact.** Building IJI per window two ways — (i) from the package's
FFT co-occurrence machinery (`.landscape_edge_counts_raster_fft` +
`.landscape_pair_count_raster_fft`) and (ii) from a brute-force moving window — agrees
to floating-point noise:

```
40x40, 4 classes, r=5:  max|FFT - brute| = 4.26e-14  (all 900 interior cells + edges)
```

This concretely confirms the whole **adjacency family** projects exactly: IJI / ent /
mutinf / relmutinf / CLUMPY are all pure functions of the same per-window co-occurrence
matrix that the shipped `contag` already FFT-projects. (Note: `fft_convolution` returns
the field transposed vs terra `[row,col]`; the shipped pipeline cancels this via
`as.vector()` + terra row-major fill at `landscape_metrics.R:1132`. A prototype that
reads the raw matrices must transpose before comparing — do not mistake this for a
projection error.)

**(B) FFT cost on small landscapes** (elapsed seconds, best of 3; composition =
`shdi`, adjacency = `contag`; absolute times vary with machine load — read the
*scaling*, not the constants):

```
 dim  nclass  radius   composition_s   adjacency_s
  50    3       5          0.03            0.07
 100    3       5          0.02            0.11
 200    3       5          0.06            0.42
  50    6       5          0.03            0.10
 100    6       5          0.01            0.28
 200    6       5          0.14            1.25
  50    3      15          0.00            0.05
 100    3      15          0.01            0.17
 200    3      15          0.09            0.39
  50    6      15          0.02            0.14
 100    6      15          0.06            0.36
 200    6      15          0.26            1.08
```

Composition is cheap (one convolution per class). Adjacency dominates because it runs
~`nclass²` pair convolutions; cost grows with cells and class count. The new adjacency
metrics add **only O(cells) arithmetic** on top of these convolutions — no extra FFTs —
so IJI / information-theory / CLUMPY are near-free once `contag`'s co-occurrence is built.

---

## 6. Implementation notes

- Landscape-level (IJI, info-theory) slot into `.msr_landscape_metrics()` /
  `.msr_adjacency_metrics()` and the `landscape_adjacency_metric_cpp` path; they are
  functions of the co-occurrence matrix already assembled there.
- Class-level (CLUMPY, PLAND/CA) need a new `class=` argument naming the focal class.
- Validate each new metric against `landscapemetrics` on circular-buffer-masked
  rasters (the approach already used for `ai`/`contag`), then check FFT projection vs
  point extraction to the documented boundary tolerance.
- Keep `relmutinf = 1` and `IJI/clumpy = NA` invariant-landscape conventions explicit
  in tests.

See also the memory note `landscape-metrics-roadmap` and the continuous-metric
counterpart `surface-metrics-roadmap`.
