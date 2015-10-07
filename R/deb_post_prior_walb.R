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

log.post.params<-function(samp, w.p, data, p, hyper, sim.data, sds){

  ## observation model
  alpha<-p["shape"]
  gamma<-p["gamma"]
  l.temp<-sim.data$L*alpha
  w.temp<-sim.data$n*gamma # shouldn't this be reproductive buffer * gamma

  llik.L<-sum(dnorm(data$L, mean=l.temp, sd=sds$L, log=TRUE))
  ##  print("here")
  llik.Negg<-sum(dnorm(data$Negg, mean=n.temp, sd=sds$Negg, log=TRUE))

  llik<-llik.L+llik.Negg

  if(length(w.p)==1) lprior<-as.numeric(log.prior.params(samp, w.p, hyper))
  else {
    lprior<-sum(log.prior.params(samp, w.p, hyper))
    ##if(!is.finite(lprior)) break
  }
  ##print(c(b, lik, prior))

  if(is.na(llik)) break
  if(is.na(lprior)) break

  return( llik + lprior )
}




log.prior.params<-function(samp, w.p, hyper){
  lp<-0
  len<-length(w.p)
  if(len==1){
    ##print(paste("w.p =", w.p, "samp =", samp[w.p[1]], sep=" "))
    lp<-list(NULL)
    names(lp)<-w.p
  }
  else{
    lp<-data.frame(matrix(0, nrow=1, ncol=len))
    names(lp)<-w.p
  }

  ##print(c(as.numeric(samp), w.p, hyper[1]))
  for(i in 1:len){
    p<-w.p[i]
    s<-as.numeric(samp[p])

    if( p == "L_m" ) "not ready yet"
    if( p == "p_Am" ) "not ready yet"
    if( p == "v" ) "not ready yet"
    if( p == "k_J" ) "not ready yet"
    if( p == "kap" ) "not ready yet"
    if( p == "T_A" ) "not ready yet"
    if( p == "T_ref" ) "not ready yet"
    if( p == "T_b" ) "not ready yet"
    if( p == "E_G" ) "not ready yet"
    if( p == "f_slope" ) lp$f_slope<-dnorm(s, mean=hyper[[p]][1], sd=hyper[[p]][2], log=TRUE)
    if( p == "f_intercept" ) lp$f_intercept<-dnorm(s, mean=hyper[[p]][1], sd=hyper[[p]][2], log=TRUE)
  }

  return(lp)

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
