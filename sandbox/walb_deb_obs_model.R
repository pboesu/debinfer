walb_obs_model_LW<-function(data, sim.data, samp){
  ec <- 1e-6
  delta_M<-samp["delta_M"]

  #length
  l.temp<-sim.data[,"L"]/delta_M

  #wet weight. this is all hard coded now, should not be!
  w_E = 23.9 # molecular weight of reserve g mol^-1
  d_v = 0.5 # specific density of structure
  mu_E = 550000 # chemical potential of reserve J / mol
  wdratio = -1.37e-3 * sim.data[,"time"] + 2.09
  omega = unname(samp['p_Am'] * w_E / (samp['v'] * d_v * mu_E))
  w.temp<-sim.data[,"L"]^3 * (1 + sim.data[,"f_n"] * omega) * d_v * wdratio # shouldn't this be reproductive buffer * gamma

  llik.L<-sum(dlnorm(data$L, meanlog=log(l.temp)+ec, sdlog=samp["sdlog.L"], log=TRUE))
  ##  print("here")
  llik.Ww<-sum(dlnorm(data$Ww, meanlog=log(w.temp)+ec, sd=samp["sdlog.Ww"], log=TRUE))

  llik<-llik.L+llik.Ww

  return(llik)
}
