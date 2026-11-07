# GDR observes two margins of one survey process. Equal removal counts or
# paired totals alone do not establish the same observed likelihood target.
.msr_gdr_selection_data <- function(model, rows) {
  .msr_gdr_native_guard(model)
  data <- model@data
  yd <- data@yDistance
  yr <- data@yRemoval
  primary <- data@numPrimary
  if (!is.matrix(yd) || !is.matrix(yr) || nrow(yd) != nrow(yr) ||
      length(primary) != 1L || !is.finite(primary) || primary < 1 ||
      primary != floor(primary) || ncol(yd) %% primary != 0L ||
      ncol(yr) %% primary != 0L || !identical(dim(data@y), dim(yr)) ||
      !isTRUE(all.equal(unname(data@y), unname(yr), tolerance = 0,
        check.attributes = FALSE)))
    stop("GDR response histories, primary layout or inherited removal mirror are inconsistent.",
      call. = FALSE)
  list(response = cbind(yd[rows, , drop = FALSE], yr[rows, , drop = FALSE]),
    protocol = list(numPrimary = primary, distance_columns = ncol(yd),
      removal_columns = ncol(yr), survey = data@survey,
      distance_breaks = data@dist.breaks, distance_units = data@unitsIn,
      removal_durations = data@period.lengths))
}
