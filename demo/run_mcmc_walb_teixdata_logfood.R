## this script runs mcmc inference on the actual data used by Teixeira et al in their paper
## run with: nohup R CMD batch run_mcmc.R &

source("R/deb_mcmc.R")
source("R/deb_solver.R")
source("R/DEB_walb_logfood.r")

## need to source this one after deb_mcmc.R to tell it to use this
## versio of the prior/posterior calcs instead of the versions in
## deb_post_prior.R
source("R/deb_post_prior_walb_logfood.R")

##library("coda")

params<-setparams.DEB.walb_logfood(f_uAsym = 1.1)
inits<-setinits.DEB.walb_logfood(from.pars = params)



#set sizestep for simulator
ss<-0.5

#load data from Teixeira paper
lit.table <- read.csv('../walb_deb/data/teixeira2014tab02.csv', skip=4)
#remove last row, because the data contains a typo
lit.table <- lit.table[-nrow(lit.table),]
#make list fit for mcmc
lit.data <-list(t = lit.table$Age.d, L = lit.table$Culmen.length.cm, Ww = lit.table$Mass.kg*1000 )

#plot data
plot(lit.data$t,lit.data$L)
plot(lit.data$t,lit.data$Ww)

#extract Tmax
Tmax <- max(lit.data$t)
abline(v=Tmax)


#solve model for published pars
teix.log.fit<-solve.DEB(walb_deb_log_food, params, inits, Tmax=Tmax, numsteps=NULL,
                    which=1, sizestep=1, verbose = FALSE)
plot(teix.log.fit)
plot(teix.log.fit, obs=data.frame(time=lit.data$t, L=lit.data$L*params['delta_M']))
## p.start is a list giving the initial guesses for the params we want
## to infer since we're not currently proposing X.h, the initial
## condition for e0 is set, so we don't need to reset the "inits"
## inits<-setinits.DEB()
hyper<-make.hypers() #make the prior slightly more vague than the std. error on the fit suggests
w.p<-c("f_uAsym", "f_lAsym", "f_rate", "f_xmid") #name the parameters that are to be estimated??

p.start<-c(1.1, 0.4, -0.03, 80) #initial values of parameters
prop.sd<-c(f_uAsym = 0.01, f_lAsym = 0.01, f_rate = 0.001, f_xmid = 5)#what is this? Metropolis-Hastings Tuning parameter?!


sds<-list(L=0.5, Ww=250)

N<-2000

log.samps<-deb.mcmc(N=N, p.start=p.start, data=lit.data, w.p=w.p, params=params, inits=inits, sim=walb_deb_log_food, sds=sds, hyper=hyper, prop.sd=prop.sd, Tmax=Tmax, cnt=50, burnin=1000, plot=TRUE, sizestep=ss, which = 1, data.times = lit.data$t, free.inits = "setinits.DEB.walb_logfood")

lit.samps<-log.samps$samps

#the output object should also store the model (name at least) and the likelihood
mcmc.out = list(samps = lit.samps, hyper = hyper, p.start = p.start, tuning = prop.sd)
save(mcmc.out, file = paste("walb_deb_logfood_samples", format(Sys.time(), format = "%Y%m%d-%H%M%S"), ".RData", sep=''))

#run multiple parallel chains
library(doParallel)
registerDoParallel(cores=4)

p.start.m <- matrix(c(1.1, 0.4, -0.03, 80,
                      1.1, 0.4, -0.03, 80,
                      1.1, 0.4, -0.03, 80,
                      1.1, 0.4, -0.03, 80),
                      ncol = 4, byrow = T)
system.time({
chains <-foreach(i=1:nrow(p.start.m)) %dopar% deb.mcmc(N=25000, p.start=p.start.m[i,], data=lit.data, w.p=w.p, params=params, inits=inits, sim=walb_deb_log_food, sds=sds, hyper=hyper, prop.sd=prop.sd, Tmax=Tmax, cnt=50, burnin=1000, plot=FALSE, sizestep=ss, which = 1, data.times = lit.data$t, free.inits = "setinits.DEB.walb_logfood")
})
save(chains, file = paste("walb_deb_logfood_samples_4chains", format(Sys.time(), format = "%Y%m%d-%H%M%S"), ".RData", sep=''))

par(mfrow=c(2,2))
plot(chains[[1]]$samps$f_uAsym, type = 'l', ylim = c(1,2))
lines(chains[[2]]$samps$f_uAsym, col = 'red')
lines(chains[[3]]$samps$f_uAsym, col = 'darkgreen')
lines(chains[[4]]$samps$f_uAsym, col = 'cornflowerblue')

plot_chains(chains, nrow = 2, ncol = 2)
pretty_pairs(do.call('rbind', (do.call('rbind', chains))))

#pdf("figs/chain_post_prior.pdf", width=12, height=8)
par(mfrow=c(1,2))
plot(lit.samps$f_lAsym,type='l', main="mcmc trace", xlab = 'iteration')

#hist(samps$f_slope, freq=FALSE)
plot(density(lit.samps$f_lAsym),xlim=c(0,2), main='f_uAsym')
lines(seq(0,2,length.out = 100),dlnorm(seq(0,2,length.out = 100), meanlog=hyper$f_lAsym[1], sdlog=hyper$f_lAsym[2]),col="red")
#legend("topleft",legend=c("prior PDF", "posterior KDE"), lty=1, col=c('red','black'))
#plot(samps$f_intercept,type='l', main="mcmc trace", xlab = 'iteration')
#plot(density(samps$f_intercept), main='f_intercept')
#lines(seq(0,4,length.out = 100),dnorm(seq(0,4,length.out = 100), mean=hyper$f_intercept[1], sd=hyper$f_intercept[2]),col="red")
#dev.off()

#solve with estimated parameters
new.params <- params
new.params['f_uAsym'] <- mean(lit.samps$f_uAsym)
new.params['f_lAsym'] <- mean(lit.samps$f_lAsym)
new.params['f_rate'] <- mean(lit.samps$f_rate)
new.params['f_xmid'] <- mean(lit.samps$f_xmid)
#new.params['f_intercept'] <- mean(samps$f_intercept[2000:5000])

new.data<-solve.DEB(walb_deb_log_food, new.params, inits = setinits.DEB.walb_logfood(from.pars = new.params), Tmax=Tmax, numsteps=NULL, which=1, sizestep=2, verbose = FALSE)

#plot(old.data)#plot debtool fit, which is the basis for the simulated data
plot(new.data)
plot(new.data, obs=data.frame(time=lit.data$t, L=lit.data$L*params['delta_M']))

get_Ww <- function(w_E = 23.9, d_v = 0.5, mu_E = 550000, sim.data=new.data, ode.pars=new.params ){
#wet weight. this is all hard coded now, should not be!
#w_E = 23.9 # molecular weight of reserve g mol^-1
#d_v = 0.5 # specific density of structure
#mu_E = 550000 # chemical potential of reserve J / mol
wdratio = -1.37e-3 * sim.data[, 'time'] + 2.09
omega = unname(ode.pars['p_Am'] * w_E / (ode.pars['v'] * d_v * mu_E))
Ww = sim.data[,'L']^3 * (1 + sim.data[,'f_n'] * omega) * d_v * wdratio
return(Ww)
}
new.Ww <- get_Ww(w_E = 23.9, d_v = 0.5, mu_E = 550000, sim.data=new.data, ode.pars=new.params )
teix.Ww <- get_Ww(sim.data = teix.fit, ode.pars=params)

#pdf("figs/post_pred.pdf", width=12, height=8)
par(mfrow=c(1,2))
plot(lit.data$t, lit.data$Ww,type='p',lwd=3,xlab = 'time (days since hatching)', ylab = 'wet weight (g)', ylim = c(0, max(lit.data$Ww)), xlim = c(0, Tmax))
lines(new.data[,"time"],new.Ww, col='red',lty=3,lwd=4)
lines(teix.fit[,"time"], teix.Ww)
legend("topleft",lwd=c(3,4,NA), lty=c(1,3,NA), pch=c(NA,NA,1), col=c('black','red','grey'),legend=c('Teixeira', 'MCMC fit','data'),bty='n')
plot(lit.data$t, lit.data$L,type='p',lwd=3,xlab = 'time (days since hatching)', ylab = 'culmen length (cm)', ylim = c(0, max(lit.data$L)), xlim = c(0, Tmax))
lines(new.data[,"time"],new.data[,"L"]/params["delta_M"], col='red',lty=3,lwd=4)
lines(teix.fit[,"time"], teix.fit[,"L"]/params["delta_M"])


#dev.off()

