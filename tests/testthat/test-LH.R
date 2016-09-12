context("Class LH")

test_that("constructor", {
  expect_is(LH(), "LH")
  expect_is(LH(method="LH axes"), "LH")
  expect_is(LH(method="LH axes", lambda=c(1.05, 1.2)), "LH")
  expect_is(pars<- sampleLH(), "data.frame")
  expect_is(LH(pars), "LH")
  
  obj<- LH()
  expect_is(LH(S3Part(obj)), "LH")
  expect_equivalent(LH(S3Part(obj)), obj)
})

test_that("subsetting", {
  obj<- LH()
  expect_is(obj[c(1,4,8),], "LH")
  expect_is(head(obj), "LH")
  expect_is(obj[1:10,], "LH")
  
  expect_is(obj$a[2], "numeric")
  expect_is(obj[[2]], "numeric")
})
