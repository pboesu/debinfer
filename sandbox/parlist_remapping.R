#untangling parlist object to map back onto old deb_mcmc
library(deBInfer)
r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                  value = 0.5, prior="norm", hypers=list(mean = 0, sd = 1),
                  prop.var=0.1, samp.type="rw")

K <- debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                  value = 5, prior="lnorm", hypers=list(meanlog = 1, sdlog = 1),
                  prop.var=0.1, samp.type="rw")

loglogsd.N <- debinfer_par(name = "loglogsd.N", var.type = "obs", fixed = FALSE,
                           value = -2, prior="norm", hypers=list(mean = 0, sd = 1),
                           prop.var=0.5, samp.type="rw")

#we also need to provide an initial condition for the DE
N <- debinfer_par(name = "loglogsd.N", var.type = "init", fixed = TRUE,
                           value = 0.1, prior="norm", hypers=list(mean = 0, sd = 1),
                           prop.var=0.5, samp.type="rw")


mcmc.pars <- setup_debinfer(r, K, loglogsd.N, N)

# do inference with deBInfer
# MCMC iterations
iter = 5000
# define burnin
burnin = 2000
# inference call

##remapping here
# get all parameter names
p.names <- sapply(mcmc.pars, function(x) x$name)
is.free <- !sapply(mcmc.pars, function(x) x$fixed)
is.init <- sapply(mcmc.pars, function(x) x$var.type)=="init"

  #inits are matched by order in deSolve --> how to ensure they are in the correct order?? add init.order variable to parameter declaration?



#get all start values
p.start <- lapply(mcmc.pars, function(x) x$value)
names(p.start) = p.names
# get the parameters that are to be estimated
w.p <- p.names[is.free]
#check what this is needed for, except for the "true" likelihood calculation
params <-

#initial values for DE (no ordering!)
inits <- sapply(mcmc.pars, function(x) x$value)[is.init]
names(inits) <- p.names[is.init]

sds <- lapply(mcmc.pars, function(x) x$value)[is.free &] # is this obsolete if it is not used in the obs model?
hyper = list
pdfs =
  prop.sd
data.times


mcmc_samples <- deb_mcmc(N=iter, p.start=p.start, data=N_obs, w.p=w.p, params=parms,
                         inits=c(N=0.1), sim=logistic_model, sds=list(N=0.01),
                         hyper=list(r=list(mean=0, sd=1),
                                    K=list(meanlog=1, sdlog=1),
                                    loglogsd.N=list(mean=-2, sd=1)),
                         pdfs = list(r='norm', K='lnorm', loglogsd.N='norm'),
                         prop.sd=c(r=0.005, K=0.1, loglogsd.N=0.5),
                         Tmax=max(N_obs$time), cnt=iter, burnin=200, plot=FALSE, sizestep=0.1, which = 1,
                         data.times = N_obs$time, obs.model=logistic_obs_model)

testthat::expect_equal(unname(colMeans(mcmc_samples$samps[burnin:iter,])/parms),c(1,1,1),tolerance = 1e-5)
