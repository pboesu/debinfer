#more utility functions for the inference outputs

#' Summary of the inference results
#'
#' A wrapper for coda::summary.mcmc
#'
#' @param object a deBInfer_result object
#' @param ... further arguments to summary.mcmc
#' @seealso \code{\link[coda]{summary.mcmc}}
#' @import coda
#' @export
summary.debinfer_result <- function(object, ...){
 summary(object$samples, ...)
}

#' is.debinfer_result
#'
#' Check debinfer_result class
#'
#' @param x an object
#' @export
is.debinfer_result <- function(x){
  if (inherits(x, "debinfer_result")) TRUE
  else FALSE
}

#' is.debinfer_parlist
#'
#' Check debinfer_parlist class
#'
#' @param x an object
#' @export
is.debinfer_parlist <- function(x){
  if (inherits(x, "debinfer_parlist")) TRUE
  else FALSE
}

#' Get starting/fixed values of DE initial values
#'
#' Accessor function for initial values
#'
#' @param x a debinfer_result or debinfer_parlist object
#' @return a named numeric vector
#' @export
deinits <- function(x){
  if (is.debinfer_result(x)){
    is.init <- vapply(x$all.params, function(x) x$var.type, character(1)) %in% c("init", "initfunc")
    inits <- vapply(x$all.params, function(x) x$value, numeric(1))[is.init]
    return(inits)
  } else {
    if (is.debinfer_parlist(x)){
      is.init <- vapply(x, function(x) x$var.type, character(1)) %in% c("init", "initfunc")
      inits <- vapply(x, function(x) x$value, numeric(1))[is.init]
      return(inits)
    } else NULL}
}

#' Get starting/fixed values of DE parameters
#'
#' Accessor function for parameters
#'
#' @param x a debinfer_result or debinfer_parlist object
#' @return a named numeric vector
#' @export
depars <- function(x){
  if (is.debinfer_result(x)){
    is.depar <- vapply(x$all.params, function(x) x$var.type, character(1))=="de"
    depars <- vapply(x$all.params, function(x) x$value, numeric(1))[is.depar]
    return(depars)
  } else {
    if (is.debinfer_parlist(x)){
      is.depar <- vapply(x, function(x) x$var.type, character(1))=="de"
      depars <- vapply(x, function(x) x$value, numeric(1))[is.depar]
      return(depars)
    } else NULL}
}

#' Get starting/fixed values of DE initial values
#'
#' Accessor function for free parameters
#'
#' @param x a debinfer_result or debinfer_parlist object
#' @return a named logical vector
#' @export
freeparams <- function(x){
  if (is.debinfer_result(x)){
    is.free <- vapply(x$all.params, function(x) !x$fixed, logical(1))
    freepars <- vapply(x$all.params, function(x) x$fixed, logical(1))[is.free]
    return(freepars)
  } else {
    if (is.debinfer_parlist(x)){
      is.free <- vapply(x, function(x) !x$fixed, logical(1))
      freepars <- vapply(x, function(x) !x$fixed, logical(1))[is.free]
      return(freepars)
    } else NULL}
}


#' Reshape posterior model solutions
#'
#' Take a list of DE model solutions and transform into a list of of matrices, one for each state variable, where each row is an iteration, and each column is a time point
#'
#' @param x a post_sim object
#' @import plyr
#' @export
reshape_post_sim <- function(x){
  if(!inherits(x, "post_sim")) stop("input not of class 'post_sim'")
  out <- list()
  out$time <- x[[1]][,'time']
  for (i in 2:ncol(x[[1]])){
    name <- colnames(x[[1]])[i]
    out[[name]] <- plyr::laply(x, function(x) x[,i])
  }
  return(out)
}
