# Coordinate mapping only. A native-model template still has to verify actual
# provider metadata and fitted/constructor blocks before calling this helper.
.msr_gdr_exact_start_map <- function(constructor_blocks, fitted_blocks, fitted_values) {
  valid <- function(x) is.character(x) && length(x) > 0L && !anyNA(x) &&
    all(grepl("^(beta|lsigma)_(lambda|phi|dist|rem|alpha|scale)$", x))
  if (!valid(constructor_blocks) || !valid(fitted_blocks))
    stop("Use explicit native fixed-regression and log-SD block identities; random modes are not constructor starts.")
  if (!is.numeric(fitted_values) || !is.null(dim(fitted_values)) ||
      length(fitted_values) != length(fitted_blocks) || any(!is.finite(fitted_values)))
    stop("Fitted values must be finite and exactly match the fitted block positions.")
  occurrence <- function(x) as.integer(stats::ave(seq_along(x), x, FUN = seq_along))
  constructor_occurrence <- occurrence(constructor_blocks)
  fitted_occurrence <- occurrence(fitted_blocks)
  key <- function(block, position) paste(block, position, sep = "#")
  target <- key(constructor_blocks, constructor_occurrence)
  original <- key(fitted_blocks, fitted_occurrence)
  if (length(target) != length(original) || !setequal(target, original))
    stop("Constructor and fitted coordinates must form the same complete block/position bijection.")
  index <- match(target, original)
  list(starts = unname(fitted_values[index]),
    map = data.frame(constructor_position = seq_along(index), block = constructor_blocks,
      within_block = constructor_occurrence, fitted_position = index),
    qualification = "coordinate_mapping_only_no_model_or_random_effect_admission")
}
