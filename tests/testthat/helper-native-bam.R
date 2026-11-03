# Shared native prediction comparison helpers. Their existing definitions
# remain unchanged; family fixtures keep separate failure-retaining caches.
bam_family_slice <- function(value, rows) {
  if (is.list(value)) { result <- lapply(value, bam_family_slice, rows = rows); attributes(result) <- attributes(value); return(result) }
  if (is.matrix(value)) result <- value[rows, , drop = FALSE] else if (length(dim(value)) == 1L) {
    labels <- dimnames(value); result <- array(as.numeric(value)[rows], length(rows),
      dimnames = if (is.null(labels)) NULL else list(labels[[1L]][rows]))
  } else result <- value[rows]
  for (name in setdiff(names(attributes(value)), c("names", "dim", "dimnames"))) attr(result, name) <- attr(value, name)
  result
}

bam_family_clean <- function(value) {
  attr(value, "prediction_contract") <- NULL
  if (is.list(value)) { result <- lapply(value, bam_family_clean); attributes(result) <- attributes(value); return(result) }
  if (is.matrix(value)) rownames(value) <- NULL else names(value) <- NULL
  value
}

bam_family_oracle <- function(entry, newdata, type, controls = list()) {
  vars <- intersect(all.vars(delete.response(terms(entry$model))), names(entry$data))
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  source <- entry$data[vars]
  result <- do.call(stats::predict, c(list(object = entry$model, newdata = rbind(source, newdata[vars]),
    type = type, se.fit = type != "lpmatrix", discrete = FALSE, n.threads = 1L), controls))
  bam_family_slice(result, nrow(source) + seq_len(nrow(newdata)))
}

