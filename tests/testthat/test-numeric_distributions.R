context("Numeric distributions")

## Check c memory management
# gctorture(on=TRUE)
# gctorture(on=FALSE)

test_that("constructor works", {
  expect_is(distri<- distriBinom(2, .6), "numericDistri")
  expect_is(distriC<- distriBinom(distri, .3), "numericDistri")
  
  res<- resC<- resS<- resP<- numeric()
  for (i in 1:1000){
    res[i]<- cumsum(distriBinom(2, .6))$cump[3]
    resP[i]<- cumsum(distri * 2)$cump[3]
    resC[i]<- cumsum(distriBinom(distri, .3))$cump[3] ## Fixed
    resS[i]<- cumsum(distri + distriC)$cump[5] ## Fixed
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
    res[i]<- cumsum(distriBinom(2, .6, log=TRUE))$cump[3]
    resP[i]<- cumsum(distri * 2)$cump[3]
    resC[i]<- cumsum(distriBinom(distri, .3, log=TRUE))$cump[3] ## Fixed
    resS[i]<- cumsum(distri + distriC)$cump[5] ## Fixed
    # print(resS[i])
  }
  
  expect_equal(unique(res), 1)
  expect_equal(unique(resP), 1)
  expect_equal(unique(resC), 1)
  expect_equal(unique(resS), 1)
  
})
