context("Class LH")

lambda<- seq(.9, 1.2, length.out=3)
broods<- 2^(0:2)
b<- 2^(0:4)
a<- seq(0.4, 0.95, length.out=3)
AFR<-  2^(0:2)

lh<- LH()

test_that("constructor", {
  expect_is(LH(), "LH")
  expect_is(LH(method="LH axes"), "LH")
  expect_is(LH(method="LH axes", lambda=c(1, 1.2)), "LH")
  expect_is(LH(method="regular"), "LH")

  expect_warning(lh<- LH(lambda=lambda, broods=broods, b=b, a=a, AFR=AFR, method="regular")) # parameters out of domain
  expect_is(lh, "LH")
  expect_equal(lh$s, lh$a) # if s is missing -> s=a
  expect_is(pars<- sampleLH(), "data.frame")
  expect_is(LH(pars), "LH")
  
  expect_is(LH(S3Part(lh)), "LH")
  expect_equivalent(LH(S3Part(lh)), lh)
})


test_that("subsetting", {
  expect_is(lh[c(1,3,5),], "LH")
  expect_is(head(lh), "LH")
  expect_is(lh[1:5,], "LH")
  
  expect_is(lh$a[2], "numeric")
  expect_is(lh[[4]], "numeric")
})


test_that("popbio", {
  expect_is(lh<- LH(method="LH axes", popbio=TRUE), "LH")
  expect_is(LH(pars=S3Part(lh), popbio=TRUE), "LH") ## test reusing popbio columns
  expect_is(LH(pars=S3Part(lh), popbio=FALSE), "LH") ## test reusing popbio columns
  
  expect_is(lh<- LH(method="regular", popbio=TRUE), "LH")
  expect_is(LH(pars=S3Part(lh), popbio=TRUE), "LH") ## test reusing popbio columns
  
  expect_warning(lh<- LH(lambda=lambda, broods=broods, b=b, a=a, AFR=AFR, method="regular", popbio=TRUE)) # parameters out of domain
  expect_is(lh, "LH")
})


test_that("plot", {
  expect_is(plot(lh), "NULL")
})
