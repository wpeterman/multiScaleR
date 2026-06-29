# profile_sigma() Held-Covariate Reuse Benchmark

Generated: 2026-06-29 09:27:32 EDT

This benchmark compares the current `profile_sigma()` implementation with a
legacy-equivalent loop that recomputes every covariate at every profile grid
point. The profile values are checked for equality before timing.

| Case | Covariates | Points per profile | Legacy median (s) | Current median (s) | Speedup |
|------|------------|--------------------|-------------------|--------------------|---------|
| `landscape_metrics_four_covariates` | 4 | 5 | 0.940 | 0.360 | 2.611x |
| `weighted_surface_metrics_four_covariates` | 4 | 5 | 1.080 | 0.450 | 2.400x |

Raw timings: `tools/benchmarks/profile_sigma_covariate_reuse_results.csv`.
