suppressPackageStartupMessages({ library(terra); library(landscapemetrics) })
source("bench_tmp/Lref.R"); options(width = 130)
mk <- function(M) { r <- terra::rast(M); terra::ext(r) <- c(0, ncol(M), 0, nrow(M)); terra::crs(r) <- "EPSG:3857"; r }
lmv <- function(f, r) { v <- suppressWarnings(tryCatch(f(r)$value, error = function(e) NA_real_)); v[1] }

scenarios <- list(
  "single class (m=1)"      = matrix(1L, 6, 6),
  "two classes (stripes)"   = matrix(rep(c(1L,2L), 18), 6, 6),
  "three classes random"    = { set.seed(5); matrix(sample(1:3, 36, TRUE), 6) },
  "NA holes (5 cls)"        = { set.seed(9); M <- matrix(sample(1:5,144,TRUE),12); M[sample(144,30)] <- NA; M },
  "rare single-cell class"  = { M <- matrix(1L, 8, 8); M[1,1] <- 2L; M[8,8] <- 3L; M }
)
metrics <- list(
  IJI       = list(ref=function(M,A,cl) ref_iji(A),               lm=lsm_l_iji),
  ent       = list(ref=function(M,A,cl) ref_info(A)["ent"],       lm=lsm_l_ent),
  mutinf    = list(ref=function(M,A,cl) ref_info(A)["mutinf"],    lm=lsm_l_mutinf),
  relmutinf = list(ref=function(M,A,cl) ref_info(A)["relmutinf"], lm=lsm_l_relmutinf),
  clumpy_c1 = list(ref=function(M,A,cl) ref_clumpy(M, cl[1]),     lm=function(r) lsm_c_clumpy(r)[1,]),
  pland_c1  = list(ref=function(M,A,cl) ref_pland(M, cl[1]),      lm=function(r) lsm_c_pland(r)[1,])
)
showna <- function(x) if (is.null(x)||length(x)==0||is.na(x)) (if (is.nan(x)) "NaN" else "NA") else sprintf("%.4f", x)
cat(sprintf("%-24s %-10s %12s %12s   %s\n", "scenario", "metric", "reference", "landscapemetrics", "match"))
for (nm in names(scenarios)) {
  M <- scenarios[[nm]]; storage.mode(M) <- "double"; r <- mk(M)
  cl <- sort(unique(stats::na.omit(as.vector(M)))); A <- adj_matrix(M, cl)
  for (mn in names(metrics)) {
    rv <- suppressWarnings(unname(metrics[[mn]]$ref(M, A, cl)))
    lv <- lmv(metrics[[mn]]$lm, r)
    eq <- (is.na(rv) && is.na(lv)) || (is.finite(rv) && is.finite(lv) && isTRUE(all.equal(rv, lv, tolerance=1e-5)))
    cat(sprintf("%-24s %-10s %12s %12s   %s\n", nm, mn, showna(rv), showna(lv), if (eq) "ok" else "DIFF"))
  }
  cat("\n")
}
cat("L2 DONE\n")
