# Independent normal-model check for the normalized tempered bridge.
# Source the larger-field runner without launching its case loop.
source("tools/benchmarks/spatial_field_cv_larger/run.R")
set.seed(644)
observed <- c(-1.2, .4, 2.1, 1.3)
prior_sd <- 1.4
noise <- .8
temperature <- seq(0, 1, by = .25)
mean_at <- function(t) t * sum(observed)/noise^2/(prior_sd^-2 + t * length(observed)/noise^2)
sd_at <- function(t) 1/sqrt(prior_sd^-2 + t * length(observed)/noise^2)
n <- 40000L
log_likelihood <- lapply(temperature, function(t) {
  states <- rnorm(n, mean_at(t), sd_at(t))
  rowSums(vapply(observed, function(y) dnorm(y, states, noise, log = TRUE), numeric(n)))
})
bridges <- lapply(seq_along(diff(temperature)), function(i)
  estimate_bridge(log_likelihood[[i]], log_likelihood[[i + 1L]], temperature[i + 1L] - temperature[i]))
bridge <- sum(vapply(bridges, function(z) z$log_ratio, 0))
V <- matrix(prior_sd^2, length(observed), length(observed))
diag(V) <- diag(V) + noise^2
L <- t(chol(V))
reference <- -(length(observed) * log(2 * pi) + sum(forwardsolve(L, observed)^2))/2 -
  sum(log(diag(L)))
stopifnot(abs(bridge - reference) < .01)
cat("Normalized Gaussian bridge difference (nats):", bridge - reference, "\n")
