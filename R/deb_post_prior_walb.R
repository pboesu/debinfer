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
  l.temp<-sim.data$l*alpha
  n.temp<-sim.data$n*gamma # shouldn't this be reproductive buffer * gamma

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

    if( p == "J.EAm" ) "not ready yet"
    if( p == "y.EX" ) "not ready yet"
    if( p == "v" ) "not ready yet"
    if( p == "k.M" ) "not ready yet"
    if( p == "g" ) lp$g<-dgamma(s, shape=hyper[[p]][1], rate=hyper[[p]][2], log=TRUE)
    if( p == "kap" ) lp$kap<-dbeta(s, hyper[[p]][1], hyper[[p]][2], log=TRUE)
    if( p == "k.J" ) "not ready yet"
    if( p == "M.PH" ) "not ready yet"
    if( p == "shape" ) "not ready yet"
    if( p == "gamma" ) "not ready yet"
    if( p == "X.h" ) "not ready yet"
    if( p == "vol" ) "not ready yet"
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
                      f_slope = NULL,
                      f_intercept = NULL){

  hyper<-list(L_m = L_m, p_Am = p_Am, v = v, k_J = k_J, kap = kap, T_A = T_A, T_ref = T_ref, T_b = T_b, E_G = E_G, f_slope = f_slope, f_intercept = f_intercept)

  return(hyper)

}
