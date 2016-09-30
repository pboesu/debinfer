# Changes in version 0.4.1.9000
##Changed functionality
* runtime plotting during the MCMC procedure `de_mcmc(... , plot=TRUE)` now omits fixed parameters
* correlation coefs in pairs plot now scale with strength of correlation 
##Minor bug fixes
* sampler messages in `de_mcmc` now follow the setting for `verbose.mcmc` rather than `verbose`
