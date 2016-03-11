## this script runs mcmc inference on the actual data used by Teixeira et al in their paper
## run with: nohup R CMD batch run_mcmc.R &

#source("R/deb_mcmc.R")
#source("R/deb_solver.R")
library("deBInfer")

#load DE model
source("demo/DEB_walb.r")

## need to source this one after deb_mcmc.R to tell it to use this
## versio of the prior/posterior calcs instead of the versions in
## deb_post_prior.R
source("demo/deb_post_prior_walb.R")

library("coda")


params<-setparams.DEB.walb()
inits<-setinits.DEB.walb(from.pars = params)


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
teix.fit<-solve_de(DEB.walb, params, inits, Tmax=Tmax, numsteps=NULL,
                    which=1, sizestep=1, verbose = FALSE)
plot(teix.fit)
plot(teix.fit, obs=data.frame(time=lit.data$t, L=lit.data$L*params['delta_M']))
## p.start is a list giving the initial guesses for the params we want
## to infer since we're not currently proposing X.h, the initial
## condition for e0 is set, so we don't need to reset the "inits"
## inits<-setinits.DEB()
#hyper<-make.hypers(kap = c(2,2), L_m = c(2.82983379, 0.08221382), p_Am = c( 6.5506723, 0.0285656), v = c(-3.50711314,  0.03332408), k_J = c(-15.02003861,   0.03332408), E_G = c(9.47104320, 0.05764439)) #make the prior slightly more vague than the std. error on the fit suggests
pdfs = list(kap = 'beta',
            L_m = 'lnorm',
            p_Am = 'lnorm',
            v = 'lnorm',
            k_J = 'lnorm',
            T_A = 'lnorm',
            T_ref = 'lnorm',
            T_b = 'lnorm',
            E_G = 'lnorm',
            f_slope = 'norm',
            f_intercept = 'norm'
            )
hyper=list(kap = list(shape1=2, shape2=2),
           L_m = list(meanlog = 2.82983379, sdlog = 0.08221382),
           p_Am = list(meanlog = 6.5506723, sdlog = 0.0285656),
           v = list(meanlog = -3.50711314,  sdlog = 0.03332408),
           k_J = list(meanlog = -15.02003861, sdlog = 0.03332408),
           E_G = list( meanlog = 9.47104320, sdlog = 0.05764439),
           f_slope = list(mean = -0.006154762,sd = 0.002985912911), #from linear fit on food data in Teixeira 2014
           f_intercept = list(mean = 1.561011905, sd = 0.4685600146)
           )
w.p<-c("f_slope", "kap", "L_m", "p_Am", "v", "k_J", "E_G", "f_intercept") #name the parameters that are to be estimated??

#to do march 6th
#write data structures for priors and hypers. done!
#in deb_mcmc split up data and prior likelihoods, so that only the former is evaluated from a user spec function


p.start<-c(-0.004, 0.9, 15, 700, 3e-2, 3e-7, 1.3e4, 1.5) #initial values of parameters; should be a named numeric, really
prop.sd<-c(f_slope=0.0001, f_intercept = 0.01, kap = 0.01, L_m = 0.2, p_Am = 10, v = 1e-3, k_J = 1e-8, E_G = 750)#what is this? Metropolis-Hastings Tuning parameter?!


sds<-list(L=0.5, Ww=250)

N<-300

lit.samps<-deb_mcmc(N=N, p.start=p.start, data=lit.data, w.p=w.p, params=params, inits=inits, sim=DEB.walb, sds=sds, hyper=hyper, pdfs = pdfs, prop.sd=prop.sd, Tmax=Tmax, cnt=30, burnin=200, plot=TRUE, sizestep=ss, which = 1, data.times = lit.data$t, free.inits = "setinits.DEB.walb")
##out<-mcmc(N=N, p.start=p.start, data, params, inits, sim=DEB1, sds, Tmax, burnin=0, cnt=50)

lit.samps<-lit.samps$samps

mcmc.out = list(samps = lit.samps, hyper = hyper, p.start = p.start, tuning = prop.sd)
save(mcmc.out, file = paste("walb_deb_teixdata_samples", format(Sys.time(), format = "%Y%m%d-%H%M%S"), ".RData", sep=''))

#loading results from file
lit.samps <- mcmc.out[['samps']]

pdf("figs/pretty_pairs_full_model_teixdat.pdf", width=12, height=8)
pretty_pairs(lit.samps[runif(10000, 2000, 30000),], trend = T)
dev.off()

#pdf("figs/chain_post_prior.pdf", width=12, height=8)
par(mfrow=c(1,2))
plot(lit.samps$f_slope,type='l', main="mcmc trace", xlab = 'iteration')
#hist(samps$f_slope, freq=FALSE)
plot(density(lit.samps$f_slope),xlim=c(-.01,0), main='f_slope')
lines(seq(-.01,0,length.out = 100),dnorm(seq(-0.01,0,length.out = 100), mean=hyper$f_slope[1], sd=hyper$f_slope[2]),col="red")
#legend("topleft",legend=c("prior PDF", "posterior KDE"), lty=1, col=c('red','black'))
#plot(samps$f_intercept,type='l', main="mcmc trace", xlab = 'iteration')
#plot(density(samps$f_intercept), main='f_intercept')
#lines(seq(0,4,length.out = 100),dnorm(seq(0,4,length.out = 100), mean=hyper$f_intercept[1], sd=hyper$f_intercept[2]),col="red")
#dev.off()

#solve with estimated parameters
new.params <- params
mean.params <- params
median.params <- params
#dut of burnin
lit.samps.nib <- lit.samps[2000:30000,]
#"f_slope", "kap", "L_m", "p_Am", "v", "k_J", "E_G"
new.params['f_slope'] <- mean(lit.samps.nib$f_slope)
for (i in w.p){
  mean.params[i] <-   mean(lit.samps.nib[,i])
  median.params[i] <-   median(lit.samps.nib[,i])
}


#new.params['f_intercept'] <- mean(samps$f_intercept[2000:5000])

new.data<-solve.DEB(DEB.walb, new.params, inits, Tmax=Tmax, numsteps=NULL,which=1, sizestep=1, verbose = FALSE)

mean.sim<-solve.DEB(DEB.walb, mean.params, inits, Tmax=Tmax, numsteps=NULL, which=1, sizestep=ss, verbose = FALSE)
median.sim<-solve.DEB(DEB.walb, median.params, inits, Tmax=Tmax, numsteps=NULL, which=1, sizestep=ss, verbose = FALSE)

#get credible intervals
#lapply



plot(mean.sim, obs=as.data.frame(median.sim))

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
new.Ww <- get_Ww(w_E = 23.9, d_v = 0.5, mu_E = 550000, sim.data=median.sim, ode.pars=median.params )
teix.Ww <- get_Ww(sim.data = teix.fit, ode.pars=params)

pdf("figs/median_pred_walb_teixdata.pdf", width=12, height=8)
par(mfrow=c(1,2))
plot(lit.data$t, lit.data$Ww,type='p',lwd=3,xlab = 'time (days since hatching)', ylab = 'wet weight (g)', ylim = c(0, max(lit.data$Ww)), xlim = c(0, Tmax))
lines(median.sim[,"time"],new.Ww, col='red',lty=2,lwd=4)
lines(teix.fit[,"time"], teix.Ww)
legend("topleft",lwd=c(3,4,NA), lty=c(1,2,NA), pch=c(NA,NA,1), col=c('black','red','grey'),legend=c('Teixeira et al. 2014', 'MCMC fit (median)','data'),bty='n')
plot(lit.data$t, lit.data$L,type='p',lwd=3,xlab = 'time (days since hatching)', ylab = 'culmen length (cm)', ylim = c(0, max(lit.data$L)), xlim = c(0, Tmax))
lines(median.sim[,"time"],median.sim[,"L"]/params["delta_M"], col='red',lty=2,lwd=4)
lines(teix.fit[,"time"], teix.fit[,"L"]/params["delta_M"])
dev.off()

par(mfrow = c(2,4))
for (i in 1:length(w.p)){
  plot(lit.samps[,i], main=w.p[i], type='l')
}


sim_from_samps <- function(sim = DEB.walb, params, samps, iteration){

  new.params <- params
  slice <- samps[iteration,]
  for (i in names(slice)) new.params[i] <- samps[iteration,i]

  new.data<-solve.DEB(sim = sim, new.params, inits, Tmax=Tmax, numsteps=NULL,
                      which=1, sizestep=ss, verbose = FALSE)

  #add step to calculate Lcul and Ww
  Ww <- get_Ww(w_E = 23.9, d_v = 0.5, mu_E = 550000, sim.data=new.data, ode.pars=new.params )
  Lcul <- new.data[,"L"]/new.params["delta_M"]

  return(cbind(new.data, Ww = Ww, Lcul = Lcul))

}

simlist <- plyr::llply(sample(nrow(lit.samps[2000:30000,]), 100), function(x)sim_from_samps(params = params, samps = lit.samps, iteration = x), .progress='text')

sima <- simplify2array(simlist)
dim(sima)

simCI <- plyr::aaply(sima, .margins = c(1,2), quantile, probs=c(0.025,0.5,0.975))

pdf("figs/simCI.pdf", width = 9, height = 5)
par(mfrow=c(2,3))
for (p in 2:ncol(simlist[[1]])){
  plot(simCI[,,3][,c(1,p)], main=dimnames(simlist[[1]])[[2]][p], type='n')

  #for (i in 1:length(simlist)){
  #  lines(simlist[[i]][,c(1,p)],lty=2,col='grey')
  #}
  for (i in 1:3){
    lines(simCI[,,i][,c(1,p)], lty=1, col=c("red"))
  }
  if (dimnames(simlist[[1]])[[2]][p] == "Ww") points(lit.data$t,lit.data$Ww)
  if (dimnames(simlist[[1]])[[2]][p] == "Lcul") points(lit.data$t, lit.data$L)
}
dev.off()
