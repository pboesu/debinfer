library(deBInfer)
context("Testing samplers")

test_that("reflecting sampler does not return out-of-bounds values", {
  for(i in rep(c(-5,-1,0,0.5,1,5),100)){
    proposal <- deBInfer:::propose_single_rev(samps=c(b=5), s.p=list(name='b', prop.var=0.7, samp.type="rw-ref01", hypers=NULL))$b
    expect_true(proposal <= 1 && proposal >= 0)
  }


})

