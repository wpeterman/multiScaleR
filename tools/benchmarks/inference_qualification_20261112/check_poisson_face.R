# Recheck profile likelihood at the reported finite theta for records whose
# independent NB2 reference could use the exact Poisson face.
source("tools/benchmarks/inference_qualification_20261112/run.R")
stopifnot(phase == "glm")
rows <- table[table$method == "profile_glm" & is.finite(table$native_gap) &
  table$native_gap > .0001, ]
results <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
  row <- rows[i, ]
  f <- simulate_case(plan[match(row$id, plan$id), ])
  initial <- native_at(f$data, f$formula, "nbinom2", FALSE, independent_nb2 = FALSE)$model
  fit <- suppressWarnings(multiScale_optim(initial, f$inputs, backend = "profile_glm",
    n_cores = 1, verbose = FALSE,
    profile_control = list(starts = matrix(c(.075, .20, .43, .71), ncol = 1L))))
  theta <- fit$dispersion_est$theta
  data <- f$data
  data$hab <- f$habitat(fit$scale_est$Mean[1L] / 1000)
  cap <- glm(f$formula, data = data, family = MASS::negative.binomial(theta))
  poisson <- glm(f$formula, data = data, family = poisson())
  data.frame(id = row$id, theta,
    finite_theta_gap = -as.numeric(logLik(fit$opt_mod)) + as.numeric(logLik(cap)),
    poisson_gap = -as.numeric(logLik(fit$opt_mod)) + as.numeric(logLik(poisson)),
    unrestricted_reference_gap = row$native_gap)
}))
utils::write.csv(results, file.path(output, "poisson-face.csv"), row.names = FALSE)
print(results, row.names = FALSE)
