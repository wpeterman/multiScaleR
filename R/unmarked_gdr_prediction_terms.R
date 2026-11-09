# GDR's native getDesign builds four fixed matrices directly. Rebind only its
# local matrix helpers, retaining its site/period/interval order and NA checks.
.msr_unmarked_gdr_frame_inputs <- function(object, type, newdata, na.rm,
    re.form, bases, raw, formulas) {
  if (!methods::is(newdata, "unmarkedFrameGDR"))
    stop("GDR frame prediction requires a native GDR frame.", call. = FALSE)
  ns <- asNamespace("unmarked")
  method <- methods::selectMethod("getDesign", class(newdata)[1L],
    fdef = methods::getGeneric("getDesign", where = ns))
  native <- .msr_unmarked_local_method(method, list())
  statements <- as.list(body(native))[-1L]
  assignment <- function(name) {
    selected <- Filter(function(x) is.call(x) && length(x) == 3L &&
      identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name(name)), statements)
    if (length(selected) != 1L) return(NULL)
    selected[[1L]][[3L]]
  }
  components <- c("lambda", "phi", "dist", "rem")
  targets <- vapply(statements, function(x) if (is.call(x) && length(x) == 3L &&
      identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]])) as.character(x[[2L]]) else "", character(1))
  count_calls <- function(x, name) {
    if (!is.call(x)) return(0L)
    as.integer(identical(x[[1L]], as.name(name))) +
      sum(vapply(as.list(x), count_calls, integer(1), name = name))
  }
  if (any(diff(match(paste0("X_", components), targets)) <= 0L, na.rm = TRUE) ||
      any(diff(match(paste0("Z_", components), targets)) <= 0L, na.rm = TRUE) ||
      count_calls(body(native), "model.frame") != 4L ||
      count_calls(body(native), "model.matrix") != 4L || count_calls(body(native), "get_Z") != 4L)
    stop("This native GDR frame evaluation order needs separate qualification.", call. = FALSE)
  fixed <- c("lam_fixed", "phi_fixed", "dist_fixed", "rem_fixed")
  data <- c("sc", "ysc", "ysc", "oc")
  for (i in seq_along(components)) {
    name <- components[i]
    expected <- paste0("model.matrix(", fixed[i], ", model.frame(", fixed[i],
      ", ", data[i], ", na.action = NULL))")
    if (!identical(assignment(fixed[i]), str2lang(paste0("reformulas::nobars(formulas$", name, ")"))) ||
        !identical(assignment(paste0("X_", name)), str2lang(expected)) ||
        !identical(assignment(paste0("Z_", name)),
          str2lang(paste0("get_Z(formulas$", name, ", ", data[i], ")"))))
      stop("This native GDR frame layout needs separate trained-design qualification.", call. = FALSE)
  }
  fixed_index <- random_index <- 0L
  fixed_frame <- function(formula, data, ...) {
    fixed_index <<- fixed_index + 1L
    if (fixed_index > 4L) stop("Unexpected native GDR fixed-design call.", call. = FALSE)
    name <- components[fixed_index]
    if (!name %in% names(bases)) return(stats::model.frame(formula, data, ...))
    structure(list(component = name, data = data), class = "msr_gdr_prediction_frame")
  }
  fixed_matrix <- function(formula, data, ...) {
    if (!inherits(data, "msr_gdr_prediction_frame")) return(stats::model.matrix(formula, data, ...))
    .msr_unmarked_fixed_design(bases[[data$component]], data$data)$X
  }
  random_matrix <- function(formula, data) {
    random_index <<- random_index + 1L
    if (random_index > 4L) stop("Unexpected native GDR random-design call.", call. = FALSE)
    name <- components[random_index]
    if (!name %in% names(bases)) return(get("get_Z", ns)(formula, data))
    .msr_unmarked_random_design(formulas[[name]], raw[[name]], data, re.form)
  }
  get_design <- .msr_unmarked_local_method(native, list(
    model.frame = fixed_frame, model.matrix = fixed_matrix, get_Z = random_matrix))
  inputs <- .msr_unmarked_local_method(get("predict_inputs_from_umf", ns),
    list(getDesign = get_design))
  value <- inputs(object, type, newdata, na.rm, re.form)
  if (fixed_index != 4L || random_index != 4L)
    stop("Native GDR component preparation did not evaluate the qualified layout.", call. = FALSE)
  value
}
