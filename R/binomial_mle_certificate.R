# Exact checks concern the represented double design, not a rounded/scaled LP.
# QR and LP results only propose witnesses. Failed exact reconstruction leaves
# a case unresolved rather than relaxing the existence criterion.
.msr_binomial_qvector <- function(x) gmp::as.bigq(as.character(x))

.msr_binomial_qinverse <- function(x) {
  q <- gmp::as.bigq(x)
  p <- nrow(x)
  R <- gmp::as.bigq(cbind(x, diag(p)))
  # Explicit exact row pivoting also handles structural zeros in factor and
  # interaction designs. No floating rank tolerance establishes invertibility.
  for (column in seq_len(p)) {
    available <- which(R[seq.int(column, p), column] != 0)
    if (!length(available)) return(NULL)
    selected <- available[[1L]] + column - 1L
    if (selected != column) {
      saved <- R[column, ]; R[column, ] <- R[selected, ]; R[selected, ] <- saved
    }
    R[column, ] <- R[column, ] / R[column, column]
    for (i in setdiff(seq_len(p), column))
      if (R[i, column] != 0) R[i, ] <- R[i, ] - R[i, column] * R[column, ]
  }
  inverse <- R[, p + seq_len(p), drop = FALSE]
  if (!all(gmp::`%*%`(q, inverse) == diag(p))) return(NULL)
  inverse
}

# Project a proposed vector into a small exact null space. Numerical active
# constraints are only proposals; callers check all original constraints.
.msr_binomial_qproject <- function(B, candidate) {
  R <- gmp::as.bigq(B); pivot <- integer(); row <- 1L
  for (column in seq_len(ncol(B))) {
    if (row > nrow(B)) break
    available <- which(R[seq.int(row, nrow(B)), column] != 0)
    if (!length(available)) next
    selected <- available[[1L]] + row - 1L
    if (selected != row) {
      saved <- R[row, ]; R[row, ] <- R[selected, ]; R[selected, ] <- saved
    }
    R[row, ] <- R[row, ] / R[row, column]
    for (i in setdiff(seq_len(nrow(B)), row))
      if (R[i, column] != 0) R[i, ] <- R[i, ] - R[i, column] * R[row, ]
    pivot <- c(pivot, column); row <- row + 1L
  }
  result <- .msr_binomial_qvector(candidate)
  free <- setdiff(seq_len(ncol(B)), pivot)
  for (i in seq_along(pivot))
    result[pivot[[i]]] <- if (length(free)) -sum(R[i, free] * result[free]) else 0
  result
}

.msr_binomial_paired_basis <- function(X, successes, trials) {
  n <- nrow(X); p <- ncol(X)
  ordering <- do.call(order, c(unname(as.list(as.data.frame(X))), list(method = "radix")))
  group <- if (n == 1L) 1L else cumsum(c(TRUE,
    rowSums(X[ordering[-1L], , drop = FALSE] != X[ordering[-n], , drop = FALSE]) != 0))
  success <- which(successes[ordering] > 0)
  failure <- which(successes[ordering] < trials[ordering])
  success <- success[!duplicated(group[success])]
  failure <- failure[!duplicated(group[failure])]
  matched <- match(group[success], group[failure])
  success <- success[!is.na(matched)]; failure <- failure[matched[!is.na(matched)]]
  if (length(success) < p) return(NULL)
  rows <- X[ordering[success], , drop = FALSE]
  scaled <- sweep(rows, 2L, pmax(apply(abs(rows), 2L, max), .Machine$double.xmin), "/")
  selected <- qr(t(scaled), LAPACK = TRUE)$pivot[seq_len(p)]
  if (is.null(.msr_binomial_qinverse(rows[selected, , drop = FALSE]))) return(NULL)
  list(kind = "paired_basis", success_rows = ordering[success[selected]],
    failure_rows = ordering[failure[selected]])
}

#' Verify a private binomial geometry certificate
#'
#' @param X,successes,trials The original represented coefficient target.
#' @param certificate A target-bound paired-row, balancing or recession witness.
#' @return A logical verification result and a reason. Missing exact arithmetic,
#'   modified inputs and malformed certificates cannot establish existence.
#' @noRd
.msr_binomial_verify_certificate <- function(X, successes, trials, certificate) {
  fail <- function(reason) list(valid = FALSE, reason = reason)
  if (!requireNamespace("gmp", quietly = TRUE)) return(fail("Optional gmp is unavailable."))
  tryCatch({
    if (!identical(certificate$design_id, digest::digest(list(X, successes, trials), algo = "sha256")))
      return(fail("The coefficient target differs from the certificate."))
    p <- ncol(X); n <- nrow(X)
    valid_rows <- function(i, size) is.numeric(i) && length(i) == size &&
      all(is.finite(i) & i == floor(i) & i >= 1 & i <= n)
    if (identical(certificate$kind, "paired_basis")) {
      s <- certificate$success_rows; f <- certificate$failure_rows
      if (!valid_rows(s, p) || !valid_rows(f, p) ||
          !all(successes[s] > 0) || !all(successes[f] < trials[f]) ||
          !all(X[s, , drop = FALSE] == X[f, , drop = FALSE]) ||
          is.null(.msr_binomial_qinverse(X[s, , drop = FALSE])))
        return(fail("Paired rows do not establish a full-rank overlap basis."))
      return(list(valid = TRUE, conclusion = "finite", arithmetic = "exact_rational"))
    }
    A <- rbind(X[successes > 0, , drop = FALSE], -X[successes < trials, , drop = FALSE])
    if (identical(certificate$kind, "positive_balance")) {
      w <- certificate$weights; rows <- certificate$corrected_rows
      if (!is.numeric(w) || length(w) != nrow(A) || any(!is.finite(w)) ||
          length(rows) != p || any(!is.finite(rows) | rows != floor(rows) | rows < 1 | rows > nrow(A)) ||
          anyDuplicated(rows) || length(certificate$corrected_weights) != p ||
          is.null(.msr_binomial_qinverse(t(A[rows, , drop = FALSE]))))
        return(fail("Balancing certificate is malformed or lacks full rank."))
      q <- gmp::as.bigq(w); q[rows] <- gmp::as.bigq(certificate$corrected_weights)
      if (!all(q > 0) || !all(gmp::crossprod(gmp::as.bigq(A), q) == 0))
        return(fail("Balancing weights are not strictly positive and exactly balanced."))
      return(list(valid = TRUE, conclusion = "finite", arithmetic = "exact_rational"))
    }
    if (identical(certificate$kind, "recession")) {
      if (!is.character(certificate$direction) || length(certificate$direction) != p)
        return(fail("Recession certificate is malformed."))
      margins <- gmp::`%*%`(gmp::as.bigq(A), gmp::as.bigq(certificate$direction))
      if (!all(margins >= 0) || !any(margins > 0))
        return(fail("The original design does not verify the recession direction."))
      return(list(valid = TRUE, conclusion = "nonexistent", arithmetic = "exact_rational",
        strict = all(margins > 0)))
    }
    fail("Unknown certificate kind.")
  }, error = function(e) fail(paste("Exact verification failed:", conditionMessage(e))))
}
