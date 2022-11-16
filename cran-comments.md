## Release summary

This is a minor update of the previous version it addresses an error raised by the vignette building process on MKL R-devel CRAN checks, for which the package was archived today (2022-11-16). My apologies for the slow response.

## Test environments
* win-builder (R-devel and R-release)
* macOS 10.13.6 High Sierra, R-release, brew (on r-hub)
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (on r-hub)
* Debian Linux, R-release, GCC (on r-hub)
* Debian Linux, R-devel, clang, ISO-8859-15 locale (on r-hub)


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Philipp H Boersch-Supan <pboesu@gmail.com>’

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1002/ece3.334
    From: man/chytrid.Rd
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1073/pnas.1311790110
    From: man/debinfer_par.Rd
    Status: 503
    Message: Service Unavailable
  URL: https://doi.org/10.1111/2041-210X.12679
    From: man/debinfer_par.Rd
          inst/CITATION
          README.md
    Status: 503
    Message: Service Unavailable

Found the following (possibly) invalid DOIs:
  DOI: 10.1111/2041-210X.12679
    From: inst/CITATION
    Status: Service Unavailable
    Message: 503
    
* All highlighted DOIs and URLs were manually checked and found to be working.
