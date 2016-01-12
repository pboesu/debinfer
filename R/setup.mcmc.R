#' setup_debinfer
#'
#' Creates an object of class debinfer_in containing initial values, parameters, prior distributions, hyperparameters
#' tuning parameters etc. to set up a debinfer analysis
#'
#' @param var.name character vector; name of the variable
#' @param var.type character vector; type of the variable "de.par" = parameter for the differential equation, "obs.par" = parameter of the observation model, "init" = initial condition for a state variable in the differential equation
#' @param fixed boolean; TRUE = parameter is taken to be fixed, FALSE = parameter is to be estimated by MCMC
#' @param value numeric; parameter value. For fixed parameters this is the value used in the analysis for free parameters this is the starting value used when setting up the MCMC chain
#' @param block integer; number of block for jont proposal; 0 or NA mean the parameter is not to be jointly proposed
#' @param distrib character; name of the probability distribution for the prior on the parameter. must conform to standard R naming i.e. beta ( foo = beta), binomial binom, Cauchy cauchy, chi-squared chisq, exponential exp, Fisher F f, gamma gamma, geometric geom, hypergeometric hyper, logistic logis, lognormal lnorm, negative binomial nbinom, normal norm, Poisson pois, Student t t, uniform unif, Weibull weibull, mvnorm
#' @param hypers list of numeric vectors, hyperparameters for the prior; mean only for mvnorm.
#' @param prop.sd numeric; tuning parameters, that is the standard deviation of the proposal distribution for each parameter
#' @param sampler character; type of sampler: "rw" random walk, "ind" = idenpendence
#'
#' @return returns an object of class debinfer_in to be fed to the mcmc function
#'
setup_debinfer <- function(var.name, var.type, fixed, value, block, distrib, hypers, prop.sd, sampler)
{

}


# > str(mcmc.p)
# List of 6
# $ :List of 5
# ..$ params: chr "fs"
# ..$ type  : chr "ind"
# ..$ var   : num 0.005
# ..$ hyp   : num [1:2] 1 1
# ..$ start : num 0.99
# $ :List of 5
# ..$ params: chr "ds"
# ..$ type  : chr "rw"
# ..$ var   : num [1:2] 3 4
# ..$ hyp   : num [1:2] 1 1
# ..$ start : num 1.71
# $ :List of 5
# ..$ params: chr "muz"
# ..$ type  : chr "rw"
# ..$ var   : num [1:2] 1 2
# ..$ hyp   : num [1:2] 5 1
# ..$ start : num 1.02
# $ :List of 5
# ..$ params: chr "eta"
# ..$ type  : chr "rw"
# ..$ var   : num [1:2] 4 5
# ..$ hyp   : num [1:2] 1 0.25
# ..$ start : num 36.4
# $ :List of 5
# ..$ params: chr "Tmin"
# ..$ type  : chr "rw"
# ..$ var   : num 0.1
# ..$ hyp   : num [1:2] 40 20
# ..$ start : num 3.9
# $ :List of 5
# ..$ params: chr "sr"
# ..$ type  : chr "ind"
# ..$ var   : num 1
# ..$ hyp   : num [1:2] 5 1
# ..$ start : num 5.11
