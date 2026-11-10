#' Native Fixed-Range GDR Response Field
#'
#' An identified extension of unmarked's GDR fit for a fixed-range exponential
#' abundance field. Use [fit_gdr_response_field()] to construct this class.
#' @slot response_field Numeric geometry, frozen component recipes, joint
#'   prediction covariance, provider/target identities and compact fit diagnostics.
#' @param object,fitted A fitted native GDR response-field object.
#' @param nsim,seed Reserved simulation controls. Field simulation remains guarded.
#' @param ... Reserved controls for guarded inherited operations.
#' @importClassesFrom unmarked unmarkedFit
#' @md
#' @exportClass unmarkedFitGDRField
#' @name unmarkedFitGDRField-class
methods::setClass("unmarkedFitGDRField", contains = "unmarkedFitGDR",
  slots = c(response_field = "list"))

