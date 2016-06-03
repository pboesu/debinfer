#more utility functions for the inference outputs

#' Summary of the inference results
#'
#' Currently just a wrapper for coda::summary.mcmc
#'
#' @param object a deBInfer_result object
#' @param ... further arguments to summary.mcmc
#' @seealso \code{\link[coda]{summary.mcmc}}
#' @import coda
#' @export
summary.debinfer_result <- function(object, ...){
 summary(object$samples, ...)
}
