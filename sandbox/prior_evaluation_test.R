#testing different ways of programatically contructing prior evaluations



get_logd_prior <- function(x, pdf, ...){
  if (pdf == 'mvnorm'){
    stop('multivariate priors not yet implemented')
  } else {
    g <- get(paste("d", pdf, sep = ""), mode = "function")
    lp <- g(x, log=TRUE, ...)
    return(lp)
  }
}


system.time(
  for (i in 1:100000) logd_prior(x=1:25, pdf='norm', hypers=c(mean=1,sd=1))
)

system.time(
  for (i in 1:100000) get_logd_prior(x=1:25, pdf='norm', mean=1, sd=1)
)

system.time(
  for (i in 1:100000) if ('norm' == 'mvnorm') { stop('multivariate priors not yet implemented')
  } else { dnorm(x=1:25, mean=1,sd=1, log=TRUE)}
)

#testing truncated prior
require(truncdist)
get_logd_prior(x=2, pdf="trunc", spec="norm", mean=1, sd=1, a=1, b=3)
#works, but might break logic?
