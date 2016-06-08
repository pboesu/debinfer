library(deSolve)
logistic_model <- function (time, y, parms) {
  with(as.list(c(y, parms)), {
    dN <- r * N * (1 - N / K)
    list(dN)
  })
}
y <- c(N = 0.1)
parms <- c(r = 0.1, K = 10)
times <- seq(0, 120, 1)
out <- ode(y, times, logistic_model, parms, method='lsoda')

## ------------------------------------------------------------------------
set.seed(143)
logistic <- as.data.frame(out[c(1,runif(35, 0, nrow(out))),]) #force include the first time-point (t=0)


## ------------------------------------------------------------------------
# add lognormal noise
parms['sdlog.N'] <- 0.05
logistic$N_noisy <- rlnorm(nrow(logistic), log(logistic$N),parms['sdlog.N'])
#observations must be ordered for solver to work
logistic <- logistic[order(logistic$time),]
names(logistic)[2]<- "N_true"

#plot(N_obs$time, N_obs$N, ylim=c(0, max(N_obs$N,N_obs$N_noisy)))
#points(N_obs$time, N_obs$N_noisy, col="red")

devtools::use_data(logistic)
