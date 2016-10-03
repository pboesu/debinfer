library(deBInfer)
context("Testing parameter constructor function")

test_that("debinfer_input checking works as expected", {
  expect_error(debinfer_par(name="p", var.type = "de", fixed = FALSE, value = 1, prior = "norm", hypers=list(mean=1, sd=1), prop.var = -0.5, samp.type = "rw"), "prop.var must be a numeric > 0")
})


test_that("prior support bound calculation works as expected", {
  param <- debinfer_par(name="p", var.type = "de", fixed = FALSE, value = 1, prior = "norm", hypers=list(mean=1, sd=1), prop.var = 0.5, samp.type = "rw")
  expect_equal(param$bounds, c(-Inf, Inf))
  param <- debinfer_par(name="p", var.type = "de", fixed = FALSE, value = 1, prior = "unif", hypers=list(min=1, max=2), prop.var = 0.5, samp.type = "rw")
  expect_equal(param$bounds, c(1, 2))
  param <- debinfer_par(name="p", var.type = "de", fixed = FALSE, value = 1, prior = "beta", hypers=list(shape1=1, shape2=1), prop.var = 0.5, samp.type = "rw")
  expect_equal(param$bounds, c(0, 1))
  param <- debinfer_par(name="p", var.type = "de", fixed = FALSE, value = 1, prior = "lnorm", hypers=list(meanlog=1, sdlog=1), prop.var = 0.5, samp.type = "rw")
  expect_equal(param$bounds, c(0, Inf))
})
