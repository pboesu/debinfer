#looking at the solver accuracy
#analytical solution
soln <- with( list(N=0.1, K=10, r=0.1),
curve(K*N/(N+(K-N)*exp(-r*t)), from = 0, to = 120, xname = "t")
)


#lsoda
points(ode(y, N_obs$time, logistic_model, parms, method='lsoda'), col="red3")
#euler on data times
lines(ode(y, N_obs$time, logistic_model, parms, method='euler'), col="green3")
#euler step size=1
lines(ode(y, 0:120, logistic_model, parms, method='euler'), col="blue3")
#euler step size 0.5
lines(ode(y, seq(0,120, by=0.5), logistic_model, parms, method='euler'), col="blue3")
#euler step size 0.1
lines(ode(y, seq(0,120, by=0.1), logistic_model, parms, method='euler'), col="blue3")
#euler step size 0.01
lines(ode(y, seq(0,120, by=0.01), logistic_model, parms, method='euler'), col="blue3")

library(microbenchmark)
microbenchmark(
  #lsoda
  ode(y, N_obs$time, logistic_model, parms, method='lsoda'),
  #euler on data times
  ode(y, N_obs$time, logistic_model, parms, method='euler'),
  #euler step size=1
  ode(y, 0:120, logistic_model, parms, method='euler'),
  #euler step size 0.5
  ode(y, seq(0,120, by=0.5), logistic_model, parms, method='euler'),
  #euler step size 0.1
  ode(y, seq(0,120, by=0.1), logistic_model, parms, method='euler'),
  #euler step size 0.01
  ode(y, seq(0,120, by=0.01), logistic_model, parms, method='euler'),
  times=10
)
