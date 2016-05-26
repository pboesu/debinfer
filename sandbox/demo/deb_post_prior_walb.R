## In this file we put the functions that calculate prior and
## posterior probabilities of the parameters. The function that
## returns the posterior probabilities must be called
## "log.post.params", and the prior functions should in
## "log.prior.params" (which will only be called from
## log.post.params). The user can update these functions so they are
## appropriate for their data and parameters. Additionally the
## function "make.hypers" makes an appropriately structured set of
## hyper parameters (which determine the prior) for the deb.mcmc
## function.



## the following functions define the posterior and prior
## distributions for the parameters of interest: J.EAm, y.EX, v, k.M,
## g, kap, k.J, M.HP, shape, gamma, X.h, vol. There is also a function
## to automatically generate the hyper parameters for the priors

log.post.params<-function(samp, w.p, data, p, pdfs, hyper, sim.data, sds, verbose.lik=FALSE){

  ##dirty fix to treat NaNs in solver output. really the model should be scaled
  ##penalty value 0 seems to lead to "wrong solutions", trying 1e99
  sim.data <- lapply(sim.data, function(x){ x[is.nan(x)] <- 1e99; return(x)})

  ode.pars <- p #variable names need to be stream lined
  ## observation model
  delta_M<-p["delta_M"]

  #length
  l.temp<-sim.data$L/delta_M

  #wet weight. this is all hard coded now, should not be!
  w_E = 23.9 # molecular weight of reserve g mol^-1
  d_v = 0.5 # specific density of structure
  mu_E = 550000 # chemical potential of reserve J / mol
  wdratio = -1.37e-3 * sim.data$time + 2.09
  omega = unname(ode.pars['p_Am'] * w_E / (ode.pars['v'] * d_v * mu_E))
  w.temp<-sim.data$L^3 * (1 + sim.data$f_n * omega) * d_v * wdratio # shouldn't this be reproductive buffer * gamma

  llik.L<-sum(dnorm(data$L, mean=l.temp, sd=sds$L, log=TRUE))
  ##  print("here")
  llik.Ww<-sum(dnorm(data$Ww, mean=w.temp, sd=sds$Ww, log=TRUE))

  llik<-llik.L+llik.Ww

  if(verbose.lik == TRUE) print(paste("lL =", llik.L, "lWw =", llik.Ww ))

  if(length(w.p)==1) lprior<-as.numeric(log_prior_params(samp, w.p, hyper))
  else {
    lprior<-sum(log_prior_params(samp, pdfs, w.p, hyper))
    ##if(!is.finite(lprior)) break
  }
  ##print(c(b, lik, prior))

  if(is.na(llik)) break
  if(is.na(lprior)) break

  return( llik + lprior )
}







## here's a function to make the hyper params in the correct format

make.hypers<-function(L_m = NULL,
                      p_Am = NULL,
                      v = NULL,
                      k_J = NULL,
                      kap = NULL,
                      T_A = NULL,
                      T_ref = NULL,
                      T_b = NULL,
                      E_G = NULL,
                      f_slope = c(-0.006154762,0.002985912911), #from linear fit on food data in Teixeira 2014
                      f_intercept = c(1.561011905,0.4685600146)){

  hyper<-list(L_m = L_m, p_Am = p_Am, v = v, k_J = k_J, kap = kap, T_A = T_A, T_ref = T_ref, T_b = T_b, E_G = E_G, f_slope = f_slope, f_intercept = f_intercept)

  return(hyper)

}
