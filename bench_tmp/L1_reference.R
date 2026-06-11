## Reference R implementations of the candidate metrics, validated against
## landscapemetrics on full small rasters. Goal: pin the exact conventions
## (log base, double-counted rook adjacency, NA handling) before FFT work.
suppressPackageStartupMessages({ library(terra); library(landscapemetrics) })
options(width = 120)

## --- double-counted rook adjacency matrix (FRAGSTATS convention) ---
adj_matrix <- function(M, classes = sort(unique(stats::na.omit(as.vector(M))))) {
  A <- matrix(0, length(classes), length(classes), dimnames = list(classes, classes))
  ix <- function(v) match(v, classes)
  acc <- function(a, b) {
    keep <- !is.na(a) & !is.na(b)
    if (!any(keep)) return(invisible())
    ia <- ix(a[keep]); ib <- ix(b[keep])
    for (j in seq_along(ia)) { A[ia[j], ib[j]] <<- A[ia[j], ib[j]] + 1; A[ib[j], ia[j]] <<- A[ib[j], ia[j]] + 1 }
  }
  if (ncol(M) > 1) acc(M[, -ncol(M), drop=FALSE], M[, -1, drop=FALSE])
  if (nrow(M) > 1) acc(M[-nrow(M), , drop=FALSE], M[-1, , drop=FALSE])
  A
}
class_counts <- function(M, classes = sort(unique(stats::na.omit(as.vector(M))))) {
  tabulate(match(M[!is.na(M)], classes), length(classes))
}

## --- reference metrics from adjacency matrix A and class counts n ---
ref_iji <- function(A) {
  m <- nrow(A); if (m < 3) return(NA_real_)
  e <- A[upper.tri(A)]; E <- sum(e); if (E == 0) return(NA_real_)
  p <- e[e > 0] / E
  (-sum(p * log(p)) / log(0.5 * m * (m - 1))) * 100
}
## information theory on normalized co-occurrence p[i,k] = A/sum(A); base = log2
ref_info <- function(A, base = 2) {
  s <- sum(A); if (s == 0) return(c(ent=NA, condent=NA, joinent=NA, mutinf=NA, relmutinf=NA))
  P <- A / s
  px <- rowSums(P)
  H <- function(v) { v <- v[v > 0]; -sum(v * log(v, base)) }
  ent <- H(px)                 # marginal entropy H(x)
  joinent <- H(as.vector(P))   # joint entropy H(x,y)
  condent <- joinent - ent     # H(y|x) (marginals equal: symmetric P)
  mutinf <- ent - condent      # I(x;y)
  relmutinf <- if (ent > 0) mutinf / ent else NA_real_
  c(ent = ent, condent = condent, joinent = joinent, mutinf = mutinf, relmutinf = relmutinf)
}
## CLUMPY (FRAGSTATS / landscapemetrics convention): pad with a -999 background
## border (full 4-neighborhood for edge cells), double-count rook adjacency,
## min_e = minimum PERIMETER (4n / 4n+2 / 4n+4).
min_perim <- function(a) { n <- trunc(sqrt(a)); mm <- a - n^2
  ifelse(mm == 0, 4*n, ifelse(a <= n*(1+n), 4*n+2, 4*n+4)) }
ref_clumpy <- function(M, focal_class) {
  P <- rbind(-999, cbind(-999, M, -999), -999)
  classes <- sort(unique(stats::na.omit(as.vector(P))))
  A <- adj_matrix(P, classes); cnt <- class_counts(P, classes)
  like <- diag(A); other <- colSums(A); mine <- min_perim(cnt)
  keep <- classes != -999
  p_i <- cnt[keep] / sum(cnt[keep])
  g_i <- like[keep] / (other[keep] - mine[keep])
  fi <- match(focal_class, classes[keep]); gi <- g_i[fi]; pi <- p_i[fi]
  if (is.nan(gi) || is.na(gi) || pi == 1) return(NA_real_)
  if (gi >= pi) (gi - pi)/(1 - pi)
  else if (pi >= 0.5) (gi - pi)/(1 - pi)
  else (gi - pi)/(-pi)
}

## ---- validate on full small rasters against landscapemetrics ----
mk_rast <- function(seed, dim = 30, nclass = 4) {
  set.seed(seed)
  r <- terra::rast(nrows = dim, ncols = dim, xmin = 0, xmax = dim, ymin = 0, ymax = dim, crs = "EPSG:3857")
  terra::values(r) <- sample(seq_len(nclass), terra::ncell(r), replace = TRUE)
  r
}
cmp <- function(seed) {
  r <- mk_rast(seed); M <- terra::as.matrix(r, wide = TRUE)
  classes <- sort(unique(as.vector(M))); A <- adj_matrix(M, classes); n <- class_counts(M, classes)
  lm_v <- function(f) f(r)$value
  info <- ref_info(A)
  data.frame(
    seed = seed,
    iji      = ref_iji(A)      - lm_v(lsm_l_iji),
    ent      = info["ent"]     - lm_v(lsm_l_ent),
    condent  = info["condent"] - lm_v(lsm_l_condent),
    joinent  = info["joinent"] - lm_v(lsm_l_joinent),
    mutinf   = info["mutinf"]  - lm_v(lsm_l_mutinf),
    relmutinf= info["relmutinf"] - lm_v(lsm_l_relmutinf),
    clumpy_c1 = ref_clumpy(M, classes[1]) - lsm_c_clumpy(r)$value[1],
    pland_c1  = (100 * n[1] / sum(n)) - lsm_c_pland(r)$value[1],
    row.names = NULL
  )
}
cat("=== reference - landscapemetrics (should be ~0) ===\n")
res <- do.call(rbind, lapply(c(1,2,3,7,11), cmp))
print(round(res, 6), row.names = FALSE)
cat("\nmax abs diff per metric:\n")
print(round(apply(abs(res[,-1]), 2, max, na.rm = TRUE), 6))
cat("L1 DONE\n")
