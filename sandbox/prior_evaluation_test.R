#testing different ways of programatically contructing prior evaluations

logd_prior <- function(x, pdf, hypers, sigma=NULL){
  if (pdf == 'mvnorm'){
    stop('multivariate priors not yet implemented')
  } else {
    lp <- do.call(paste("d",pdf, sep=''), args=append(list(x, log=TRUE), hypers))
    return(lp)
  }
}

get_logd_prior <- function(x, pdf, ...){ #ellipsis and get more efficient, but how to pass in named hypers in the place of ellipsis?
  if (pdf == 'mvnorm'){
    stop('multivariate priors not yet implemented')
  } else {
    g <- get(paste("d", pdf, sep = ""), mode = "function")
    lp <- g(x, log=TRUE, ...)
    return(lp)
  }
}


match_logd_prior <- function(x, pdf, ...){ #ellipsis and get more efficient, but how to pass in named hypers in the place of ellipsis?
  if (pdf == 'mvnorm'){
    stop('multivariate priors not yet implemented')
  } else {
    g <- match.fun(paste("d", pdf, sep = ""))
    lp <- g(x, log=TRUE, ...)
    return(lp)
  }
}
##what about assembling call from pdf & hypers at sart of mcmc, then evaluating call?


microbenchmark::microbenchmark(
  do.call = logd_prior(x=1, pdf='norm', hypers=c(mean=1,sd=1)),
  get = get_logd_prior(x=1, pdf='norm', mean=1, sd=1),
  match.fun = match_logd_prior(x=1, pdf='norm', mean=1, sd=1),
  times=1e5
)


#testing truncated prior
require(truncdist)
get_logd_prior(x=2, pdf="trunc", spec="norm", mean=1, sd=1, a=1, b=3)

logd_prior(x=2, pdf="trunc2", hypers = list(spec="norm", mean=1, sd=1, a=1, b=3))
#works, but might break logic?

dtrunc2 <- truncdist::dtrunc

