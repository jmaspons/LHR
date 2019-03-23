context("Class Numeric distributions (S3)")

## Check c memory management
# gctorture(on=TRUE)
# gctorture(on=FALSE)

test_that("constructor", {
  expect_is(distri<- distriBinom(2, .6), "numericDistri")
  expect_is(distriC<- distriBinom(distri, .3), "numericDistri")
  
  res<- resC<- resS<- resP<- numeric()
  for (i in 1:1000){
    res[i]<- cumP(distriBinom(2, .6))$cump[3]
    resP[i]<- cumP(distri * 2)$cump[3]
    resC[i]<- cumP(distriBinom(distri, .3))$cump[3] ## Fixed
    resS[i]<- cumP(distri + distriC)$cump[5] ## Fixed
    # print(resS[i])
  }
  
  expect_equal(unique(res), 1)
  expect_equal(unique(resP), 1)
  expect_equal(unique(resC), 1)
  expect_equal(unique(resS), 1)
  
  ## logP
  distri<- distriBinom(2, .6, log=TRUE)
  distriC<- distriBinom(distri, .3, log=TRUE)
  
  res<- resC<- resS<- resP<- numeric()
  for (i in 1:1000){
    res[i]<- cumP(distriBinom(2, .6, log=TRUE))$cump[3]
    resP[i]<- cumP(distri * 2)$cump[3]
    resC[i]<- cumP(distriBinom(distri, .3, log=TRUE))$cump[3] ## Fixed
    resS[i]<- cumP(distri + distriC)$cump[5] ## Fixed
    # print(resS[i])
  }
  
  expect_equal(unique(res), 1)
  expect_equal(unique(resP), 1)
  expect_equal(unique(resC), 1)
  expect_equal(unique(resS), 1)
  
})

test_that("methods", {
  distri<- distriBinom(2, .6)
  distriC<- distriBinom(distri, .3)
  
  expect_is(mean(distri), "numeric")
  expect_is(var(distri), "numeric")
  expect_is(summary(distri), "data.frame")
  expect_is(quantile(distri), c("numeric", "integer"))
  expect_is(sdistri(distri), "data.frame")
  expect_is(ddistri(0:2, distri), "numeric")
  expect_is(pdistri(0:2, distri), "numeric")
  expect_is(qdistri(seq(0, 1, length=5), distri), c("numeric", "integer"))
  expect_is(rdistri(10, distri), c("numeric", "integer"))
})

