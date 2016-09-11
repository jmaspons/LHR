context("Class discretePopSim (S3)")

test_that("discretePopSim generics", {
  pop<- discretePopSim(broods=1, b=1, j=.5, a=.5, breedFail=0,
                       varJ=0, varBreedFail=0, seasonVar=1,
                       sexRatio=NA, matingSystem=NA, N0=2, replicates=15, tf=10)

  expect_is(summary(pop, dt = 1), "data.frame")
  expect_is(r(pop, dt = 1), "matrix")
  expect_is(lambda(pop, dt = 1), "matrix")
  expect_is(trendsProp(pop, dt = 1), "data.frame")
  
  expect_is(G(pop), "numeric")
  expect_is(Gmean(pop), "numeric")
})
