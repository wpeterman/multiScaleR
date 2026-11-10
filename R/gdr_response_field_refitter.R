# Preserve the ordinary callback interface while keeping the field target and
# trained component recipes distinct from the iid native constructor.
.msr_gdr_field_refitter <- function(model, cold_starts, audit_dir, provenance,
    tmb_settings) {
  field <- .msr_gdr_field_guard(model)
  settings <- field$settings
  requested <- list(method = settings$method, se = settings$se,
    threads = settings$threads, control = settings$dots$control)
  requested[names(settings$dots)[names(settings$dots) != "control"]] <-
    settings$dots[names(settings$dots) != "control"]
  requested[names(settings$policy$requested)] <- settings$policy$requested
  if (!is.list(tmb_settings) || !setequal(names(tmb_settings), names(requested)) ||
      !identical(tmb_settings[names(requested)], requested))
    stop("Field refits must explicitly retain the original qualified tmb_settings.", call. = FALSE)
  plain <- function(x) {
    if (is.list(x) && !is.object(x)) return(all(vapply(x, plain, logical(1))))
    is.null(x) || (is.atomic(x) && !is.object(x) && is.null(dim(x)) && length(x) <= 64L)
  }
  if (!is.list(provenance) || !plain(provenance) ||
      length(serialize(provenance, NULL, version = 2L)) > 65536L)
    stop("Use only compact plain provenance for field refits.", call. = FALSE)
  original_frame <- model@data
  formulas <- lapply(model@formlist, deparse)
  original_key <- .msr_gdr_hash(list(field$prediction_sha256, model@opt$par,
    field$field_id, formulas, model@K, model@mixture, model@keyfun))
  force(field); force(cold_starts); force(audit_dir); force(original_frame)
  callback <- function(model, data, context = list(return_model = TRUE)) {
    current <- .msr_gdr_field_guard(model)
    key <- .msr_gdr_hash(list(current$prediction_sha256, model@opt$par,
      current$field_id, lapply(model@formlist, deparse), model@K, model@mixture, model@keyfun))
    if (!identical(key, original_key) || !methods::is(data, "unmarkedFrameGDR"))
      stop("The GDR field refitter belongs to another original fitted context.", call. = FALSE)
    covariates <- context$covs
    if (is.null(covariates)) covariates <- character()
    if (!is.character(covariates) || anyNA(covariates) || anyDuplicated(covariates) ||
        any(!covariates %in% names(original_frame@siteCovs)))
      stop("Field refits require explicit original static site covariates.", call. = FALSE)
    .msr_gdr_spatial_guard(model, covariates)
    retained <- data; retained@siteCovs <- original_frame@siteCovs
    if (!identical(retained, original_frame) ||
        !identical(names(data@siteCovs), names(original_frame@siteCovs)) ||
        !identical(rownames(data@siteCovs), rownames(original_frame@siteCovs)) ||
        !identical(data@siteCovs[setdiff(names(data@siteCovs), covariates)],
          original_frame@siteCovs[setdiff(names(original_frame@siteCovs), covariates)]))
      stop("A GDR field refit changed histories, survey protocol, row order or nonspatial columns.", call. = FALSE)
    result <- .msr_gdr_field_fit(model, data, cold_starts, audit_dir, field,
      .objective_only = isFALSE(context$return_model), .provenance = provenance)
    result@call <- model@call
    result@response_field$provenance <- provenance
    result
  }
  scope <- list2env(list(field = field, cold_starts = cold_starts, audit_dir = audit_dir,
    original_frame = original_frame, original_key = original_key, provenance = provenance),
    parent = environment(.msr_gdr_field_refitter))
  environment(callback) <- scope
  callback
}
