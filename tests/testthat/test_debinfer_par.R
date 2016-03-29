library(deBInfer)
context("Testing parameter constructor function")

test_that("debinfer_input checking works as expected", {
  expect_error(debinfer_par(name="p", var.type = "de", fixed = FALSE, value = 1, prior = "norm", hypers=list(mean=1, sd=1), prop.var = -0.5, samp.type = "rw"), "prop.var must be a numeric > 0")
})
