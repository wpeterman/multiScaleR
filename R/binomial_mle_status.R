#' Check private binomial coefficient existence
#'
#' @param X Finite numeric coefficient design, including any intercept column.
#' @param successes Integer success counts in the design's row order.
#' @param trials Positive integer trials for each row. Binary data use one.
#' @param tolerance Positive tolerance for proposing numerical witnesses. Exact
#'   verification of the original design determines acceptance independently.
#' @param .solver Internal LP provider used only for isolated failure tests.
#' @return A conservative status with a target-bound exact certificate when
#'   finite existence or nonexistence is established. No status grants selection
#'   eligibility. Missing optional arithmetic and unresolved cases stay unknown.
#' @noRd
.msr_binomial_mle_status <- function(X, successes, trials = rep(1, nrow(X)),
    tolerance = 1e-8, .solver = NULL) {
  if (!is.matrix(X) || !is.numeric(X) || !nrow(X) || !ncol(X) ||
      any(!is.finite(X)) || !is.numeric(successes) || length(successes) != nrow(X) ||
      !is.numeric(trials) || length(trials) != nrow(X) ||
      any(!is.finite(successes) | !is.finite(trials) | trials < 1 |
        trials != floor(trials) | successes < 0 | successes > trials | successes != floor(successes)) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0)
    stop("Finite-MLE checks require a finite design and aligned integer binomial counts.")
  provider_version <- if (is.null(.solver)) "not_attempted" else "injected_test_provider"
  design_id <- digest::digest(list(X, successes, trials), algo = "sha256")
  out <- function(status, reason, ...) list(status = status, reason = reason,
    qualification = "exact_original_design_certificate", solver = "lpSolve",
    solver_version = provider_version, design_id = design_id,
    tolerance = tolerance, selection_eligible = FALSE, ...)
  A0 <- rbind(X[successes > 0, , drop = FALSE], -X[successes < trials, , drop = FALSE])
  column_scale <- apply(abs(A0), 2L, max)
  if (nrow(A0) < ncol(A0) || any(column_scale == 0) || anyDuplicated(A0, MARGIN = 2L))
    return(out("non_identifiable", "The coefficient design is rank deficient."))
  if (!requireNamespace("gmp", quietly = TRUE))
    return(out("unknown", "Optional gmp exact arithmetic is unavailable."))
  bind <- function(certificate) { certificate$design_id <- design_id; certificate }
  verify <- function(certificate) .msr_binomial_verify_certificate(X, successes, trials, certificate)
  # p independent rows supporting both outcomes rule out every nonzero
  # recession direction without an LP or expansion into individual trials.
  paired <- .msr_binomial_paired_basis(X, successes, trials)
  if (!is.null(paired)) {
    paired <- bind(paired)
    if (isTRUE(verify(paired)$valid))
      return(out("finite", "An exact full-rank paired-row certificate establishes overlap.",
        certificate = paired, witness = list(kind = "paired_basis")))
  }
  A <- sweep(A0, 2L, column_scale, "/")
  singular <- svd(A, nu = 0L, nv = 0L)$d
  lower <- min(singular) - 64 * .Machine$double.eps * max(dim(A)) * norm(A, "F")
  if (!is.finite(lower) || lower <= 0)
    return(out("unknown", "A stable witness basis could not be proposed.", rank_status = "unresolved"))
  if (is.null(.solver)) {
    if (!requireNamespace("lpSolve", quietly = TRUE))
      return(out("unknown", "Optional lpSolve is unavailable."))
    provider_version <- as.character(utils::packageVersion("lpSolve")); .solver <- lpSolve::lp
  }
  if (!is.function(.solver)) stop("The LP provider must be a function.")
  solve <- function(...) tryCatch(.solver(...), error = function(e)
    list(status = NA_integer_, error = conditionMessage(e)))
  n <- nrow(A); p <- ncol(A)
  valid <- function(z, size) is.list(z) && is.numeric(z$status) && length(z$status) == 1L &&
    is.finite(z$status) && z$status == 0 &&
    is.numeric(z$solution) && length(z$solution) == size && all(is.finite(z$solution))
  dual <- solve(direction = "max", objective.in = c(rep(0, n), 1),
    const.mat = rbind(cbind(t(A), colSums(A)), c(rep(1, n), n)),
    const.dir = rep("=", p + 1L), const.rhs = c(rep(0, p), 1), scale = 0)
  w <- NULL
  if (valid(dual, n + 1L)) {
    w <- dual$solution[seq_len(n)] + dual$solution[n + 1L]
    if (all(w > 0)) {
      # Use column-pivoted QR to propose a minor, not the leading columns of
      # a wide default QR. Exact inversion still establishes its actual rank.
      rows <- qr(t(A * w), LAPACK = TRUE)$pivot[seq_len(p)]
      inverse <- .msr_binomial_qinverse(t(A0[rows, , drop = FALSE]))
      if (!is.null(inverse)) {
        q <- gmp::as.bigq(w)
        correction <- gmp::`%*%`(inverse, gmp::crossprod(gmp::as.bigq(A0), q))
        corrected <- q[rows] - .msr_binomial_qvector(correction)
        certificate <- bind(list(kind = "positive_balance", weights = w,
          corrected_rows = rows, corrected_weights = as.character(corrected)))
        if (isTRUE(verify(certificate)$valid)) {
          q[rows] <- corrected; weights <- as.numeric(q)
          return(out("finite", "Strictly positive weights balance the original design exactly.",
            certificate = certificate, witness = list(weights = weights,
              minimum_weight = min(weights), positivity_margin = min(weights),
              exact_strict_positivity = TRUE), dual_status = dual$status))
        }
      }
    }
  }
  primal <- solve(direction = "max", objective.in = c(colSums(A), -colSums(A)),
    const.mat = rbind(cbind(A, -A), rep(1, 2L * p)),
    const.dir = c(rep(">=", n), "<="), const.rhs = c(rep(0, n), 1), scale = 0)
  direction_certificate <- function(z) {
    if (!valid(z, 2L * p) && !valid(z, 2L * p + 1L)) return(NULL)
    scaled <- gmp::as.bigq(z$solution[seq_len(p)]) - gmp::as.bigq(z$solution[p + seq_len(p)])
    direction <- scaled / gmp::as.bigq(column_scale)
    certificate <- bind(list(kind = "recession", direction = as.character(direction)))
    if (isTRUE(verify(certificate)$valid)) return(certificate)
    near <- which(abs(drop(A %*% as.numeric(scaled))) <= tolerance)
    if (length(near)) {
      basis <- qr(t(A[near, , drop = FALSE]))
      near <- near[basis$pivot[seq_len(min(basis$rank, p))]]
      if (length(near)) {
        direction <- .msr_binomial_qproject(A0[near, , drop = FALSE], direction)
        certificate$direction <- as.character(direction)
        if (isTRUE(verify(certificate)$valid)) return(certificate)
      }
    }
    NULL
  }
  certificate <- direction_certificate(primal); strict <- NULL
  if (!is.null(certificate)) {
    check <- verify(certificate)
    if (!check$strict) {
      # Maximizing the sum can choose a non-strict direction even when a strict
      # separator exists. Only a separately verified strict witness proves that subtype.
      strict <- solve(direction = "max", objective.in = c(rep(0, 2L * p), 1),
        const.mat = rbind(cbind(A, -A, -1), c(rep(1, 2L * p), 0)),
        const.dir = c(rep(">=", n), "<="), const.rhs = c(rep(0, n), 1), scale = 0)
      proposed <- direction_certificate(strict)
      if (!is.null(proposed) && isTRUE(verify(proposed)$strict)) {
        certificate <- proposed; check <- verify(certificate)
      }
    }
    separation <- if (check$strict) "complete" else "unresolved_subtype"
    if (!check$strict && !is.null(w)) {
      support <- which(w > 0)
      if (length(support) && length(support) <= p + 1L) {
        balance <- .msr_binomial_qproject(t(A0[support, , drop = FALSE]), gmp::as.bigq(w[support]))
        if (all(balance >= 0) && any(balance > 0) &&
            all(gmp::crossprod(gmp::as.bigq(A0[support, , drop = FALSE]), balance) == 0))
          separation <- "quasi_complete"
      }
    }
    return(out("nonexistent", "An exact recession direction rules out a finite coefficient maximum.",
      separation = separation, certificate = certificate,
      witness = list(direction = stats::setNames(as.numeric(gmp::as.bigq(certificate$direction)), colnames(X))),
      dual_status = dual$status, primal_status = primal$status,
      solver_errors = Filter(Negate(is.null), list(dual$error, primal$error, strict$error))))
  }
  out("unknown", "No exact witness was established; convergence cannot certify existence.",
    dual_status = dual$status, primal_status = primal$status,
    solver_errors = Filter(Negate(is.null), list(dual$error, primal$error)))
}
