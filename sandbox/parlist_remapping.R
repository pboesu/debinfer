#untangling parlist object to map back onto old deb_mcmc

#testthat::skip_on_cran()
# define logistic model
logistic_model <- function(time, y, parms) {
  with(as.list(c(y, parms)), {
    dN <- r * N * (1 - N / K)
    list(dN)
  })
}
expect_true(exists("logistic_model"))
# set initial value for simulation
y <- c(N = 0.1)
# set parameter values
parms <- c(r = 0.1, K = 10)
# set simulation time points
times <- seq(0, 120, 1)
# solve ODE
out <- ode(y, times, logistic_model, parms, method='lsoda')
# sample from simulated data
set.seed(143)
N_obs <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) #force include the first time-point (t=0)
# add lognormal noise
parms['loglogsd.N'] <- -4.6
N_obs$N_noisy <- rlnorm(nrow(N_obs), log(N_obs$N),exp(parms['loglogsd.N']))
# observations must be ordered for solver to work
N_obs <- N_obs[order(N_obs$time),]

# define an observation model
# NB: lognormal errors are not great really for the obs model - should be changed to sth that actually allows N to be zero instead of using epsilon correction
logistic_obs_model<-function(data, sim.data, sds, samp){

  llik.N<-sum(dlnorm(data$N_noisy, meanlog=log(sim.data$N+1e-6), sdlog=exp(samp[['loglogsd.N']]), log=TRUE))

  llik<-llik.N

  return(llik)
}

library(deBInfer)
r <- debinfer_par(name = "r", var.type = "de", fixed = FALSE,
                  value = 0.5, prior="norm", hypers=list(mean = 0, sd = 1),
                  prop.var=0.1, samp.type="rw")

K <- debinfer_par(name = "K", var.type = "de", fixed = FALSE,
                  value = 5, prior="lnorm", hypers=list(meanlog = 1, sdlog = 1),
                  prop.var=0.1, samp.type="rw")

loglogsd.N <- debinfer_par(name = "loglogsd.N", var.type = "obs", fixed = FALSE,
                           value = -10, prior="norm", hypers=list(mean = 0, sd = 1),
                           prop.var=0.5, samp.type="rw")

#we also need to provide an initial condition for the DE
N <- debinfer_par(name = "N", var.type = "init", fixed = TRUE,
                           value = 0.1)


mcmc.pars <- setup_debinfer(r, K, loglogsd.N, N)

# do inference with deBInfer
# MCMC iterations
iter = 5000
# define burnin
burnin = 2000
# inference call

samps <- de_mcmc(N = iter, data=N_obs, de.model=logistic_model, obs.model=logistic_obs_model, all.params=mcmc.pars,
                    Tmax = max(N_obs$time), data.times=N_obs$time, cnt=50, burnin=0.1,
                    plot=TRUE, sizestep=0.1, which=1)

samps <- de_mcmc(N = 10000, data=N_obs, de.model=logistic_model, obs.model=logistic_obs_model, all.params=mcmc.pars,
                 Tmax = max(N_obs$time), data.times=N_obs$time, cnt=50, burnin=0.1,
                 plot=TRUE, sizestep=0.1, w.t=1, which=2,
                 myswitch=NULL,
                 mymap=NULL, verbose =FALSE)







############ manual remapping

##remapping here
# get all parameter names
p.names <- sapply(mcmc.pars, function(x) x$name)
is.free <- !sapply(mcmc.pars, function(x) x$fixed)
is.init <- sapply(mcmc.pars, function(x) x$var.type)=="init"

  #inits are matched by order in deSolve --> how to ensure they are in the correct order?? add init.order variable to parameter declaration?



#get all start values
p.start <- lapply(mcmc.pars, function(x) x$value)[is.free]
names(p.start) = p.names[is.free]
# get the parameters that are to be estimated
w.p <- p.names[is.free]
#check what this is needed for, except for the "true" likelihood calculation
params <- unlist(lapply(mcmc.pars, function(x) x$value))
names(params) <-  p.names
#initial values for DE (no ordering!)
inits <- sapply(mcmc.pars, function(x) x$value)[is.init]
names(inits) <- p.names[is.init]

sds <- NULL# is this obsolete if it is not used in the obs model?
hyper = lapply(mcmc.pars, function(x) x$hyper)
names(hyper) <- p.names
pdfs = lapply(mcmc.pars, function(x) x$prior)
names(pdfs) = p.names

prop.sd <- sapply(mcmc.pars, function(x) x$prop.var)[is.free]
names(prop.sd) <- p.names[is.free]
#hardcoded/userinput
data.times = N_obs$time
Tmax=max(N_obs$time)


mcmc_samples <- deb_mcmc(N=iter, p.start=p.start, data=N_obs, w.p=w.p, params=parms,
                         inits=inits, sim=logistic_model, sds=sds,
                         hyper=hyper,
                         pdfs = pdfs,
                         prop.sd=prop.sd,
                         Tmax=Tmax, cnt=1, burnin=200, plot=TRUE, sizestep=0.1, which = 1,
                         data.times = N_obs$time, obs.model=logistic_obs_model, verbose = TRUE)

testthat::expect_equal(unname(colMeans(mcmc_samples$samps[burnin:iter,])/parms),c(1,1,1),tolerance = 1e-5)
