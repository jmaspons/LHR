context("Class LH")

lambda<- seq(.9, 1.2, by=0.1)
broods<- 2^(0:2)
b<- 2^(0:4)
a<- seq(0.4, 0.95, length.out=5)
AFR<- 1:5

test_that("constructor", {
  expect_is(LH(), "LH")
  expect_is(LH(method="LH axes"), "LH")
  expect_is(LH(method="LH axes", lambda=c(1, 1.2)), "LH")
  expect_is(LH(method="regular"), "LH")

  expect_is(LH(lambda=lambda, broods=broods, b=b, a=a, AFR=AFR, method="regular"), "LH")
  
  expect_is(pars<- sampleLH(), "data.frame")
  expect_is(LH(pars), "LH")
  
  obj<- LH()
  expect_is(LH(S3Part(obj)), "LH")
  expect_equivalent(LH(S3Part(obj)), obj)
})


test_that("subsetting", {
  obj<- LH()
  expect_is(obj[c(1,3,5),], "LH")
  expect_is(head(obj), "LH")
  expect_is(obj[1:5,], "LH")
  
  expect_is(obj$a[2], "numeric")
  expect_is(obj[[4]], "numeric")
})


test_that("popbio", {
  expect_is(LH(method="LH axes", popbio=TRUE), "LH")
  expect_is(LH(method="regular", popbio=TRUE), "LH")
  expect_is(LH(lambda=lambda, broods=broods, b=b, a=a, AFR=AFR, method="regular", popbio=TRUE), "LH")
})
