OS X builds of R packages for R-devel fail because of missing gfortran-4.8

I am trying to build an R package (github.com/pboesu/debinfer) that has dependencies that require the compilation of Fortran code. `r: release` builds work fine with both linux and OSX, but builds with `r: devel` fail with OSX because the Fortran code in the dependencies cant be compiled:

https://travis-ci.org/pboesu/debinfer/jobs/158581638#L594
```
gfortran-4.8   -fPIC  -g -O2  -c daux.f -o daux.o
make: gfortran-4.8: No such file or directory
make: *** [daux.o] Error 1
ERROR: compilation failed for package ‘deSolve’
* removing ‘/Users/travis/R/Library/deSolve’
```
This happens both for the default OSX image (https://travis-ci.org/pboesu/debinfer/jobs/158495404#L585) and the `xcode7.3` image (above). Eventhough the latter is supposedly the one that resembles the CRAN check setup (according to https://cran.r-project.org/web/checks/check_flavors.html)

I'd be grateful for any pointers on how to resolve this issue.
