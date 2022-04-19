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
{
  parlist <- list(...)
  if (!all(sapply(parlist, class) %in% c("debinfer_par", "debinfer_cov")))
    stop("input arguments need to be of class debinfer_par")
  #check for joint proposals
  #for each unique cov matrix, check that dimensions and dimension names match names and number of associated parameters
  names(parlist) <- vapply(parlist, function(x)
    x$name, character(1))
  structure(parlist, class = "debinfer_parlist")
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
      if (pdf == 'trunc'){#handle inconsistent behaviour of truncdist, which silently fails (i.e. produces incorrect output) with log = TRUE
        lp <- log(do.call(paste("d",pdf, sep=''), args=append(list(x, log = FALSE), hypers)))
      } else {
        lp <- do.call(paste("d",pdf, sep=''), args=append(list(x, log=TRUE), hypers))
      }
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
#' @param prior character; name of the probability distribution for the prior on the parameter. Must conform to standard R naming of d/r function pairs, e.g. beta ( foo = beta), binomial binom, Cauchy cauchy, chi-squared chisq, exponential exp, Fisher F f, gamma gamma, geometric geom, hypergeometric hyper, logistic logis, lognormal lnorm, negative binomial nbinom, normal norm, Poisson pois, Student t t, uniform unif, Weibull weibull. Priors from the truncdist package are available by default. User priors can be provided but must be available in the environment from which de_mcmc is called.
#' @param hypers list of numeric vectors, hyperparameters for the prior; mean only for mvnorm. Can include trunc for truncated pdfs from package truncdist.
#' @param prop.var numeric; tuning parameters. For Normal proposals (`samp.type="rw"` or `samp.type="rw-ref"`), this must be a positive number representing the standard deviation of the proposal distribution for each parameter. For the asymmetric uniform proposal distribution (`samp.type="rw-unif"`) two positive numeric values are required and the proposal will then have the bounds `prop.var[1]/prop.var[2]*current_proposal` and `prop.var[2]/prop.var[1]*current_proposal`. See Boersch-Supan et al. (2016).
#' @param samp.type character; type of sampler: "rw" = Normal random walk, "ind" = independence, "rw-unif" = asymmetric uniform distribution, "rw-ref" = reflecting random walk sampler on the bounds of the prior support (cf. Hoff 2009, Chapter 10.5.1; Yang and Rodriguez 2013)
#'
#' @return returns an object of class debinfer_par to be fed to the mcmc setup function
#' @references Boersch-Supan et al. 2016, MEE 8:511-518 \url{https://doi.org/10.1111/2041-210X.12679}
#'
#'             Hoff 2009, A First Course in Bayesian Statistical Methods, Springer
#'
#'             Yang and Rodriguez 2013, PNAS 110:19307-19312 \url{https://doi.org/10.1073/pnas.1311790110}
#' @export
debinfer_par <- function(name, var.type, fixed, value, joint=NULL, prior=NULL, hypers=NULL, prop.var=NULL, samp.type=NULL){
  #check inputs
  if(!is.character(name)) stop("name must be of type character")
  if(!var.type %in% c("de","obs","init", "initfunc")) stop('var.type must be one of c("de","obs","init", "initfunc")')
  if(var.type != "initfunc"){
    if(!is.logical(fixed)) stop("fixed must be boolean")
    if(!is.numeric(value)) stop("value must be numeric")
    if(!fixed & (is.null(prior) | is.null(hypers) | is.null(prop.var) | is.null(samp.type))) stop("free parameters require a specification of prior, hypers, prop.var and samp.type")
    if(fixed & !(is.null(prior) | is.null(hypers) | is.null(prop.var) | is.null(samp.type))) warning(paste(name, "is treated as a fixed parameters. Ignoring prior, hypers, prop.var and samp.type specification."))
    if(!fixed) if(!samp.type %in% c("rw", "rw-unif","ind", "rw-ref")) stop('samp.type must be one of c("rw", "rw-unif","ind", "rw-ref")')
    if(!fixed) if(samp.type == "rw") if(!is.numeric(prop.var) | prop.var < 0 | length(prop.var)!=1) stop("prop.var must be a numeric > 0 of length 1 for sampler type 'rw'")
    if(!fixed) if(samp.type == "rw-ref") if(!is.numeric(prop.var) | prop.var < 0 | length(prop.var)!=1) stop("prop.var must be a numeric > 0 of length 1 for sampler type 'rw-ref'")
    if(!fixed) if(samp.type == "rw-ref") if(prop.var >= 1) warning("prop.var should be << 1 for efficient sampling with sampler type 'rw-ref'")
    if(!fixed) if(samp.type == "rw-unif") if(!is.numeric(prop.var) | all(prop.var < 0) | length(prop.var)!=2) stop("prop.var must be a numeric > 0 of length 2 for sampler type 'rw-unif'")
    if(!fixed) if(samp.type == "rw-unif") if(prop.var[1] >= prop.var[2])stop("prop.var[1] must be smaller than prop.var[2] for sampler type 'rw-unif'")
    #checks for prior and hypers?
    if(!is.null(joint)) if(!is.character(joint)) stop("joint needs to be of type character (name of covariance matrix)")
    #get limits of prior support for reflection sampler
    if(!fixed){
      bounds <-  do.call(paste("q", prior, sep=""), c(list(p = c(0,1)), hypers))
    } else {
      bounds <- NA
    }
  } else {
    #check initfunc
  }


  structure(list(name = name,
                 var.type = var.type,
                 fixed = fixed,
                 value = value,
                 joint = joint,
                 prior = prior,
                 bounds = bounds,
                 hypers = hypers,
                 prop.var = prop.var,
                 samp.type = samp.type), class = "debinfer_par")
}

#' debinfer_cov
#'
#' @param var.names names of the parameters that are to be proposed together
#' @param sigma covariance matrix
#' @param samp.type character; type of sampler. currently only "rw" = Normal random walk is implemented for multivariate proposals
#' @param name name of the joint block
#'
#' @return a debinfer_cov object
#' @export
#'
debinfer_cov <- function(var.names, sigma=diag(length(names)), name , samp.type = "rw"){
  if(!is.character(var.names)) stop("var.names must be a character vector")
  if(!is.matrix(sigma) | !is.numeric(sigma)) stop("sigma must be a numeric matrix")
  if (!isSymmetric(sigma)) stop("sigma must be symmetric")
  if (!all(eigen(sigma, only.values = TRUE)$values>0)) warning("sigma does not appear to be positive semi-definite")
  if(any(dim(sigma)!=length(var.names))) stop("length(var.names) does not match dimensions of sigma")
  colnames(sigma)<-var.names
  rownames(sigma)<-var.names
  structure(list(sigma=sigma, name = name, samp.type = samp.type), class="debinfer_cov")
}

