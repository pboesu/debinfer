# maybe this function should create two objects, one with the free parameters, that is then used as a template for the samples array, and one with the fixed values, that is only evaluated during the make.states call and the likelihood calculations

#' setup_debinfer
#'
#' Creates an object of class debinfer_parlist containing initial values, parameters, prior distributions, hyperparameters
#' tuning parameters etc. to set up a debinfer analysis
#'
#' @param ... debinfer_par objects to be combined into a debinfer_parlist
#' @return returns an S3 object of class debinfer_parlist to be fed to the mcmc function
#' @export
setup_debinfer <- function(...)
{ parlist <- list(...)
  if (!all(sapply(parlist, class) == "debinfer_par")) stop("input arguments need to be of class debinfer_par")
  names(parlist)<-sapply(parlist, function(x)x$name)
  structure(parlist, class="debinfer_parlist")
}




#' logd_prior
#'
#' Evaluates the log probability density of value given a name of a prior pdf and the corresponding hyperparameters
#'
#' @param x numeric; vector of values.
#' @param pdf character; name of a probability function. Must conform to base R nomenclature of d/r function pairs. Can include trunc for truncated pdfs from package truncdist.
#' @param hypers list; a list of parameters to be passed to the density function.
#'
#' @return the value of the log density function evaluated at \code{x}
#' @importFrom truncdist dtrunc
#' @export
#'
#'
logd_prior <- function(x, pdf, hypers){
  if (pdf == 'mvnorm'){
    stop('multivariate priors not implemented')
    } else {
    lp <- do.call(paste("d",pdf, sep=''), args=append(list(x, log=TRUE), hypers))
    return(lp)
    }
}


#' debinfer_par
#'
#' Creates an object containing all the necessary bits for a parameter i.e. initial values, prior distributions,
#' hyper-parameters, tuning parameters, etc. to set up a debinfer analysis
#'
#' @param name character vector; name of the variable
#' @param var.type character vector; type of the variable "de" = parameter for the differential equation, "obs" = parameter of the observation model, "init" = initial condition for a state variable in the differential equation
#' @param fixed boolean; TRUE = parameter is taken to be fixed, FALSE = parameter is to be estimated by MCMC
#' @param value numeric; parameter value. For fixed parameters this is the value used in the analysis for free parameters this is the starting value used when setting up the MCMC chain
#' @param joint integer; number of block for joint proposal; NULL means the parameter is not to be jointly proposed
#' @param prior character; name of the probability distribution for the prior on the parameter. must conform to standard R naming of d/r function pairs, i.e. beta ( foo = beta), binomial binom, Cauchy cauchy, chi-squared chisq, exponential exp, Fisher F f, gamma gamma, geometric geom, hypergeometric hyper, logistic logis, lognormal lnorm, negative binomial nbinom, normal norm, Poisson pois, Student t t, uniform unif, Weibull weibull, mvnorm
#' @param hypers list of numeric vectors, hyperparameters for the prior; mean only for mvnorm. Can include trunc for truncated pdfs from package truncdist.
#' @param prop.var numeric; tuning parameters, that is the standard deviation of the proposal distribution for each parameter
#' @param samp.type character; type of sampler: "rw" = Normal random walk, "ind" = independence, "rw-unif" = asymmetric uniform distribution
#'
#' @return returns an object of class debinfer_par to be fed to the mcmc setup function
#' @export
debinfer_par <- function(name, var.type, fixed, value, joint=NULL, prior=NULL, hypers=NULL, prop.var=NULL, samp.type=NULL){
  #check inputs
  if(!is.character(name)) stop("name must be of type character")
  if(!var.type %in% c("de","obs","init")) stop('var.type must be one of c("de","obs","init")')
  if(!is.logical(fixed)) stop("fixed must be boolean")
  if(!is.numeric(value)) stop("value must be numeric")
  if(!fixed & (is.null(prior) | is.null(hypers) | is.null(prop.var) | is.null(samp.type))) stop("free parameters require a specification of prior, hypers, prop.var and samp.type")
  if(fixed & !(is.null(prior) | is.null(hypers) | is.null(prop.var) | is.null(samp.type))) warning(paste(name, "is treated as a fixed parameters. Ignoring prior, hypers, prop.var and samp.type specification."))
  if(!fixed) if(!samp.type %in% c("rw", "rw-unif","ind")) stop('samp.type must be one of c("rw", "rw-unif","ind)')
  if(!fixed) if(samp.type == "rw") if(!is.numeric(prop.var) | prop.var < 0 | length(prop.var)!=1) stop("prop.var must be a numeric > 0 of length 1 for sampler type 'rw'")
  if(!fixed) if(samp.type == "rw-unif") if(!is.numeric(prop.var) | all(prop.var < 0) | length(prop.var)!=2) stop("prop.var must be a numeric > 0 of length 2 for sampler type 'rw-unif'")
  if(!fixed) if(samp.type == "rw-unif") if(prop.var[1] >= prop.var[2])stop("prop.var[1] must be smaller than prop.var[2] for sampler type 'rw-unif'")
  #checks for prior and hypers?
  if(!is.null(joint)) stop("joint proposals are not yet implemented")

  structure(list(name = name,
              var.type = var.type,
              fixed = fixed,
              value = value,
              joint = joint,
              prior = prior,
              hypers = hypers,
              prop.var = prop.var,
              samp.type = samp.type), class = "debinfer_par")
}
