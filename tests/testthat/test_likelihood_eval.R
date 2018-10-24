library(deBInfer)
context("Testing functions used in the likelihood evaluation for prior and posterior")

test_that("logd_prior calculates correct log densities", {
  expect_equal(logd_prior(1, 'norm', hypers=list(mean=2,sd=1)), dnorm(1, mean=2,sd=1,log=TRUE))
  expect_equal(logd_prior(1, 'norm', hypers=list(mean=2,sd=1)), log(dnorm(1, mean=2,sd=1,log=FALSE)))
  expect_equal(logd_prior(1:10, 'norm', hypers=list(mean=2,sd=1)), dnorm(1:10, mean=2,sd=1,log=TRUE))
  expect_equal(logd_prior(1:10, 'gamma', hypers=list(shape=2,scale=1)), dgamma(1:10, shape=2,scale=1,log=TRUE))
})

test_that("logd_prior works with truncdist", {
  expect_equal(logd_prior(1:10, 'trunc', hypers=list(spec='norm', mean=2,sd=1, a=0.5)), log(dtrunc(1:10, spec='norm', mean=2, sd=1, a=0.5)))
  expect_equal(logd_prior(0:10, 'trunc', hypers=list(spec='norm', mean=2,sd=1, a=0.5)), log(dtrunc(0:10, spec='norm', mean=2, sd=1, a=0.5)))
  #test untruncated dtrunc returns same result as stats equivalent
  expect_equal(logd_prior(0:10, 'trunc', hypers=list(spec='norm', mean=2,sd=1)), dnorm(0:10, mean=2, sd=1, log = TRUE))
})
