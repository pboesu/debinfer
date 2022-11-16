# Changes in version 0.4.4
##Minor bug fixes and improvements
* `pairs.debinfer_result` now handles chains stuck at the initial value more gracefully.
# Changes in version 0.4.3
##Minor bug fixes and improvements
* clarifications in documentation
* patched bug in `truncdist::dtrunc`
* internal type checking and control flow tidied up to be compatible with upcoming changes in R 4.2.0
# Changes in version 0.4.2
##New Functionality
* New sampler type `rw-ref`, an implementation of the reflecting random walk described in Yang and Rodriguez 2013. The sampler reflects on the bounds of the associated prior distribution.
##Changed functionality
* runtime plotting during the MCMC procedure `de_mcmc(... , plot=TRUE)` now omits fixed parameters
* correlation coefs in pairs plot now scale with strength of correlation 
##Minor bug fixes
* sampler messages in `de_mcmc` now follow the setting for `verbose.mcmc` rather than `verbose`
* added github URLs to DESCRIPTION
* updated CITATION to reflect publication of the deBInfer application paper
* internal function log_post_params now gracefully handels NaNs in the likelihood 
