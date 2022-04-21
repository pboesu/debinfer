#to avoid PDF compression warnings from cran checks, build pkg on command line with compression flag
# R CMD build --compact-vignettes="gs+qpdf" debinfer

#then do cran like checks from source tarball via r-hub with particulare checks for 4.2.0 changes to length 1 vs length > 1 in logical tests
rhub::check(path = '../deBInfer_0.4.3.tar.gz', platform = c("debian-clang-devel", "windows-x86_64-devel", "macos-highsierra-release-cran", "ubuntu-gcc-devel"), env_vars = c(`_R_CHECK_LENGTH_1_LOGIC2_` = "verbose"), check_args = '--as-cran')
