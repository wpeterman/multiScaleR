# Surface roughness prefix-cache assessment

The optional cache calculates unweighted surface roughness by accumulating
sample variance in distance order. The ordinary calculation recomputes the
sample SD for each requested radius. Roundoff can affect which point a
continuous optimizer reports on a flat section of the likelihood.

## Design and results

The paired runs used 24 datasets with 48 buffers each. Eight independent seeds
were crossed with ordinary elevations, 15% missing cells, and elevations
shifted by 10 billion units. Distances included ties, and cells on a radius
were included. The same native Gaussian model was fitted with the cache off
and on using both the continuous search and a declared 13-radius grid. Run
order alternated within search and seed. An independent base R sample SD was
also evaluated at each grid radius.

All 96 fits completed. All metric missing-value masks agreed. The maximum
absolute difference from base R was 1.77e-12. All 24 grid fits selected the
same radius with matching candidate status, with a maximum candidate
log-likelihood difference of 1.60e-11. Four of the eight large-datum
continuous searches returned slightly different radii, with a largest
difference of 2.94e-10 map units. Across all fits, the largest absolute
difference in fitted response was 8.34e-13.

| Search | Median scalar seconds | Median cached seconds |
|:--|--:|--:|
| Continuous | 0.990 | 0.925 |
| Declared grid | 0.645 | 0.500 |

The timings describe small buffers during a shared Windows workload. They
support keeping the prefix algorithm as an opt-in implementation for further
performance study. The default scalar path is unchanged. This study covers one
roughness predictor and one response model. Search equivalence across other
likelihoods, multiple roughness predictors, and nearly tied grid maxima
requires separate evidence.

## Reproduce

Run from the source root with the development DLL already built.

```text
Rscript tools/benchmarks/surface_sq_cache_20261112/run.R OUTPUT
```

[Metric results](metrics.csv), [paired fit results](fits.csv), and
[session information](session.txt) are included.
