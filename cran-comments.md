## Release summary

This is a resubmission. Tests and vignettes have been further optimised. Check times on win-builder are now below 10 minutes (518 seconds for R-release, and 489 seconds for R-devel).

Otherwise this remains a minor update of the previous version it addresses warnings raised by R-devel CRAN checks about the use of && and || with vectors of length >1, and type checking with class(x), rather than inherits().

## Test environments
* local macOS 10.14.6, R 4.1.3
* win-builder (R-devel and R-release)
* Ubuntu Linux 20.04.1 LTS, R-release, gcc (on r-hub)
* Debian linux , R-devel (2022-04-20 r82221), clang (on r-hub)


## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Philipp H Boersch-Supan <pboesu@gmail.com>'

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
