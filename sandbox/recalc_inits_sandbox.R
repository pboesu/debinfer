#sandbox file to test updating of initial values
library(deBInfer)
#Lotka-Volterra model from http://strimas.com/r/lotka-volterra/
# parameters
pars <- c(alpha = 1, beta = 0.2, delta = 0.5, gammap = 0.2)
# initial state
init <- c(Nx = 1, Ny = 2)
# times
times <- seq(0, 100, by = 1)

lv_model <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    d_x <- alpha * Nx - beta * Nx * Ny
    d_y <- delta * beta * Nx * Ny - gammap * Ny
    return(list(c(Nx = d_x, Ny = d_y)))
  })
}
lv_results <- ode(init, times, lv_model, pars)
plot(lv_results)

#sample data
data_times <- seq(0,100, by = 3)
lv_data <- subset(as.data.frame(lv_results), time %in% data_times)

pars["sdlog.N"] <- 0.05
lv_data$x_noisy <- rlnorm(nrow(lv_data), log(lv_data$Nx),pars["sdlog.N"])
lv_data$y_noisy <- rlnorm(nrow(lv_data), log(lv_data$Ny),pars["sdlog.N"])

plot(lv_data$time, lv_data$y, type = "o", ylim = c(0,max(lv_data$Nx, lv_data$Ny)))
points(lv_data$time, lv_data$y_noisy, pch = "+")
lines(lv_data$time, lv_data$x, type = "o", col = "red")
points(lv_data$time, lv_data$y_noisy, pch = "+", col = "red")

# the observation model
lv_obs_model <- function(data, sim.data, samp){
  llik.x <- sum(dlnorm(data$x_noisy, meanlog = log(sim.data[,"Nx"] + 1e-6),
                       sdlog = samp[["sdlog.N"]], log = TRUE))
  llik.y <- sum(dlnorm(data$y_noisy, meanlog = log(sim.data[,"Ny"] + 1e-6),
                       sdlog = samp[["sdlog.N"]], log = TRUE))
  llik.N = llik.x + llik.y
  return(llik.N)
}

#regular inference
alph <- debinfer_par(name = "alpha", var.type = "de", fixed = FALSE,
                  value = 0.9, prior = "norm", hypers = list(mean = 0, sd = 1),
                  prop.var = 0.00001, samp.type="rw")
bet <- debinfer_par(name = "beta", var.type = "de", fixed = FALSE,
                  value = 0.21, prior = "norm", hypers = list(mean = 0, sd = 1),
                  prop.var = 0.00001, samp.type = "rw")
delt <- debinfer_par(name = "delta", var.type = "de", fixed = TRUE,
                    value = pars['delta'])
gam <- debinfer_par(name = "gammap", var.type = "de", fixed = TRUE,
                    value = pars['gammap'])
sdlog.N <- debinfer_par(name = "sdlog.N", var.type = "obs", fixed = FALSE,
                        value = 0.05, prior = "lnorm", hypers = list(meanlog = 0, sdlog = 1),
                        prop.var = c(3,4), samp.type = "rw-unif")
Nx <- debinfer_par(name = "Nx", var.type = "init", fixed = FALSE, value = 0.5, prior = "lnorm", hypers = list(meanlog = 0, sdlog = 1),
                   prop.var = c(4,5), samp.type = "rw-unif")
Ny <- debinfer_par(name = "Ny", var.type = "init", fixed = TRUE, value = 2)

mcmc.pars <- setup_debinfer(alph, bet, delt, gam, sdlog.N, Nx, Ny)

#run inference
iter <- 50
mcmc_samples <- de_mcmc(N = iter, data = lv_data, de.model = lv_model,
                        obs.model = lv_obs_model, all.params = mcmc.pars,
                        Tmax = max(lv_data$time), data.times = lv_data$time, cnt = 500,
                        plot = FALSE, verbose.mcmc = TRUE, solver = "ode")

plot(mcmc_samples)
#post_prior_densplot(mcmc_samples)


#post_traj <- post_sim(mcmc_samples, n = 200, times = seq(0,100, by = 0.1), burnin = 200, output = 'all', prob = 0.95)
#plot(post_traj, plot.type = 'ensemble', col = "#FF000040")

#now let's redefine Ny
recalc_Ny <- function(inits, params){
  inits["Ny"] <- inits["Nx"] + 1
  return(inits)
}

Ny <- debinfer_par(name = "Ny", var.type = "initfunc", value = NA, fixed = TRUE, deinitfunc = recalc_Ny)

mcmc.pars2 <- setup_debinfer(alph, bet, delt, gam, sdlog.N, Nx, Ny)


iter <- 2500
mcmc_samples2 <- de_mcmc(N = iter, data = lv_data, de.model = lv_model,
                        obs.model = lv_obs_model, all.params = mcmc.pars2,
                        Tmax = max(lv_data$time), data.times = lv_data$time, cnt = 500,
                        plot = FALSE, verbose.mcmc = TRUE, solver = "ode")

plot(mcmc_samples2)
post_traj2 <- post_sim(mcmc_samples2, n = 200, times = seq(0,100, by = 0.1), burnin = 200, output = 'all', prob = 0.95)
plot(post_traj2, plot.type = 'ensemble', col = "#FF000040")
points(lv_data$time, lv_data$y_noisy, pch = "+", col = "black")
