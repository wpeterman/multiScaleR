## De-risk the CLUMPY closed form before implementing:
##   like_i  = diag of double-counted co-occurrence (= 2 x physical i-i edges)
##   other_i = 4 * n_i              (every cell side is like / unlike / edge-background)
##   g_i     = like_i / (other_i - min_perim_i)
## Claim: this equals the validated ref_clumpy (which uses colSums on a -999-padded
## matrix) AND landscapemetrics, on complete rasters and on circular-buffer windows.
suppressPackageStartupMessages({ library(terra); library(landscapemetrics) })
source("bench_tmp/Lref.R"); options(width = 120)

## closed-form clumpy from a matrix M (NA = outside window / background)
closed_clumpy <- function(M, focal_class) {
  classes <- sort(unique(stats::na.omit(as.vector(M))))
  A <- adj_matrix(M, classes)                 # in-window double-counted co-occurrence
  n <- class_counts(M, classes)
  fi <- match(focal_class, classes); if (is.na(fi)) return(NA_real_)
  like_i  <- unname(diag(A)[fi])
  other_i <- 4 * n[fi]                         # <-- the closed form
  min_i   <- min_perim(n[fi])
  g_i <- like_i / (other_i - min_i)
  p_i <- n[fi] / sum(n)
  if (p_i == 1) return(NA_real_)
  if (g_i >= p_i) (g_i - p_i)/(1 - p_i)
  else if (p_i >= 0.5) (g_i - p_i)/(1 - p_i)
  else (g_i - p_i)/(-p_i)
}

mk <- function(M) { r <- terra::rast(M); terra::ext(r) <- c(0, ncol(M), 0, nrow(M)); terra::crs(r) <- "EPSG:3857"; r }

## (1) complete rasters vs landscapemetrics, all classes
cat("=== (1) complete rasters: closed-form vs landscapemetrics (per class) ===\n")
for (seed in c(1, 4, 9)) {
  set.seed(seed); M <- matrix(sample(1:4, 900, TRUE), 30); storage.mode(M) <- "double"
  r <- mk(M); lm <- lsm_c_clumpy(r)
  for (k in sort(unique(as.vector(M)))) {
    cf <- closed_clumpy(M, k); lv <- lm$value[lm$class == k]
    cat(sprintf("  seed=%d class=%d  closed=% .6f  lm=% .6f  %s\n",
                seed, k, cf, lv, if (isTRUE(all.equal(cf, lv, tol = 1e-6))) "ok" else "DIFF"))
  }
}

## (2) circular-buffer windows (the multiScaleR context): the validation target is
##     landscapemetrics run ON the masked circular window (same approach as ai/contag).
##     Compare BOTH closed-form (other=4n) and ref_clumpy (boundary-skipped) to lm.
cat("\n=== (2) circular-buffer windows vs lsm_c_clumpy on the masked window ===\n")
cat(sprintf("%-16s %5s %12s %12s %12s\n", "window", "class", "closed(4n)", "ref(skip)", "lm_window"))
set.seed(7); big <- matrix(sample(1:4, 60*60, TRUE), 60); radius <- 8
ok_closed <- ok_ref <- nwin <- 0
for (ctr in list(c(20,20), c(30,42), c(15,33), c(45,25))) {
  rr <- ctr[1]; cc <- ctr[2]
  ri <- (rr-radius):(rr+radius); ci <- (cc-radius):(cc+radius)
  sub <- big[ri, ci]; dmat <- sqrt(outer((ri-rr)^2, (ci-cc)^2, "+")); sub[dmat > radius] <- NA
  storage.mode(sub) <- "double"; rw <- mk(sub); lmw <- suppressWarnings(lsm_c_clumpy(rw))
  for (k in sort(unique(stats::na.omit(as.vector(sub))))) {
    cf <- closed_clumpy(sub, k); rf <- ref_clumpy(sub, k)
    lv <- lmw$value[lmw$class == k]; if (length(lv) == 0) lv <- NA_real_
    nwin <- nwin + 1
    ok_closed <- ok_closed + isTRUE(all.equal(cf, lv, tol = 1e-6))
    ok_ref    <- ok_ref + isTRUE(all.equal(rf, lv, tol = 1e-6))
    cat(sprintf("(%2d,%2d)%10s %5d % 12.6f % 12.6f % 12.6f\n", rr, cc, "", k, cf, rf, lv))
  }
}
cat(sprintf("\nmatches to lm-on-window:  closed(4n) = %d/%d   ref(skip) = %d/%d\n",
            ok_closed, nwin, ok_ref, nwin))
cat("L4 DONE\n")
