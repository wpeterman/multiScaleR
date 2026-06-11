## Validated reference implementations (match landscapemetrics 2.2.1 exactly).
adj_matrix <- function(M, classes = sort(unique(stats::na.omit(as.vector(M))))) {
  A <- matrix(0, length(classes), length(classes), dimnames = list(classes, classes))
  ix <- function(v) match(v, classes)
  acc <- function(a, b) { keep <- !is.na(a) & !is.na(b); if (!any(keep)) return(invisible())
    ia <- ix(a[keep]); ib <- ix(b[keep])
    for (j in seq_along(ia)) { A[ia[j], ib[j]] <<- A[ia[j], ib[j]] + 1; A[ib[j], ia[j]] <<- A[ib[j], ia[j]] + 1 } }
  if (ncol(M) > 1) acc(M[, -ncol(M), drop=FALSE], M[, -1, drop=FALSE])
  if (nrow(M) > 1) acc(M[-nrow(M), , drop=FALSE], M[-1, , drop=FALSE])
  A
}
class_counts <- function(M, classes = sort(unique(stats::na.omit(as.vector(M))))) {
  tabulate(match(M[!is.na(M)], classes), length(classes))
}
ref_iji <- function(A) { m <- nrow(A); if (m < 3) return(NA_real_)
  e <- A[upper.tri(A)]; E <- sum(e); if (E == 0) return(NA_real_)
  p <- e[e > 0] / E; (-sum(p * log(p)) / log(0.5 * m * (m - 1))) * 100 }
ref_info <- function(A, base = 2) { s <- sum(A)
  if (s == 0) return(c(ent=NA, condent=NA, joinent=NA, mutinf=NA, relmutinf=NA))
  P <- A / s; px <- rowSums(P); H <- function(v) { v <- v[v > 0]; -sum(v * log(v, base)) }
  ent <- H(px); joinent <- H(as.vector(P)); condent <- joinent - ent
  mutinf <- ent - condent
  # landscapemetrics convention: relmutinf = 1 when mutual information is 0
  relmutinf <- if (mutinf == 0) 1 else mutinf / ent
  c(ent = ent, condent = condent, joinent = joinent, mutinf = mutinf, relmutinf = relmutinf) }
min_perim <- function(a) { n <- trunc(sqrt(a)); mm <- a - n^2
  ifelse(mm == 0, 4*n, ifelse(a <= n*(1+n), 4*n+2, 4*n+4)) }
ref_clumpy <- function(M, focal_class) {
  P <- rbind(-999, cbind(-999, M, -999), -999)
  classes <- sort(unique(stats::na.omit(as.vector(P))))
  A <- adj_matrix(P, classes); cnt <- class_counts(P, classes)
  like <- diag(A); other <- colSums(A); mine <- min_perim(cnt); keep <- classes != -999
  p_i <- cnt[keep] / sum(cnt[keep]); g_i <- like[keep] / (other[keep] - mine[keep])
  fi <- match(focal_class, classes[keep]); if (is.na(fi)) return(NA_real_)
  gi <- g_i[fi]; pi <- p_i[fi]
  if (is.nan(gi) || is.na(gi) || pi == 1) return(NA_real_)
  if (gi >= pi) (gi - pi)/(1 - pi) else if (pi >= 0.5) (gi - pi)/(1 - pi) else (gi - pi)/(-pi) }
ref_pland <- function(M, focal_class) { n <- class_counts(M); classes <- sort(unique(stats::na.omit(as.vector(M))))
  fi <- match(focal_class, classes); if (is.na(fi)) return(0); 100 * n[fi] / sum(n) }
