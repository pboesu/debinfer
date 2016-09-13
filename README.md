[![Travis-CI Build Status](https://travis-ci.org/pboesu/debinfer.svg?branch=master)](https://travis-ci.org/pboesu/debinfer)
[![codecov](https://codecov.io/gh/pboesu/debinfer/branch/master/graph/badge.svg)](https://codecov.io/gh/pboesu/debinfer)
[![DOI](https://zenodo.org/badge/21877/pboesu/debinfer.svg)](https://zenodo.org/badge/latestdoi/21877/pboesu/debinfer)


## deBInfer: Bayesian inference for dynamical models of biological systems in R

1. Differential equations (DEs) are commonly used to model the temporal evolution of biological systems, but statistical    methods for comparing DE models to data and for parameter inference are relatively poorly developed.  This is especially problematic in the context of biological systems where observations are often noisy and only a small number of time points may be available.
2. Bayesian approaches offer a coherent framework for parameter inference that can account for multiple sources of uncertainty, while making use of prior information. We present deBInfer, an R package implementing a Bayesian framework for parameter inference in DEs. This approach offers a rigorous methodology for parameter inference as well as modeling the link between unobservable model states and parameters, and observable quantities. 
3. deBInfer  provides templates for the DE model, the observation model and data likelihood, and the model parameters and their prior distributions. A Markov chain Monte Carlo (MCMC) procedure processes these inputs to estimate the posterior distributions of the parameters and any derived quantities, including the model trajectories. Further functionality is provided to facilitate MCMC diagnostics and the visualisation of the posterior distributions of model parameters and trajectories. 
4.  The templating approach makes deBInfer applicable to a wide range of DE models and we demonstrate its application to ordinary and delay DE models for population ecology. 

For more information read our [preprint](https://arxiv.org/abs/1605.00021) or get in touch with pboesu@gmail.com

Software development is supported by [NSF grant PLR-1341649](http://www.nsf.gov/awardsearch/showAward?AWD_ID=1341649).
