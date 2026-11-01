path_component_fixture <- function(kernel, family, response, support, design) {
  f <- private_sq_fixture(kernel, family, response = response)
  if (support == "unique") {
    d <- f$inputs$d_list[[2L]]; nearest <- min(d)
    d[2L] <- (nearest + min(d[d > nearest])) / 2
    f$inputs$d_list[[2L]] <- d
  }
  f$data$band <- factor(as.integer(f$data$site) <= 2)
  f$scales[[1L]] <- switch(design, canonical = ~1 + (1 | year),
    contrast = ~1 + temp + (1 | year), cell_means = ~0 + band + (1 | year))
  .msr_scale_prepare(f$formula, f$inputs, f$data, f$scales, family, "key", reference_scale = 250)
}

