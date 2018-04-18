## Release summary

This is a minor update of the previous version. Most importantly it addresses a warning raised by R-devel CRAN checks about break being used in the wrong context.

## Test environments
* local Windows install R 3.4.4
* win-builder (R-devel and R-release)
* ubuntu 14.04.5 (on travis-ci), R 3.4.4 and R-devel (2018-04-18 r74614)
* OS X 10.12.6 (on travis-ci) R 3.4.4



## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE on win-builder R-release:

* New submission

  Possibly mis-spelled words in DESCRIPTION:
  unobservable (15:31)

* The note highlighted a possibly mis-spelled word, which is in fact a correctly spelled English adjective.
