context("Class discretePopSim (S3)")

test_that("discretePopSim generics", {
  pop<- discretePopSim(broods=1, b=1, j=.5, a=.5, breedFail=0,
                       varJ=0, varBreedFail=0, seasonVar=1,
                       sexRatio=NA, matingSystem=NA, N0=2, replicates=15, tf=10)

  expect_is(summary(pop, dt = 1), "numeric")
  expect_is(r(pop, dt = 1), "matrix")
  expect_is(lambda(pop, dt = 1), "matrix")
  expect_is(trendsProp(pop, dt = 1), "numeric")
  
  expect_is(G(pop), "numeric")
  expect_is(Gmean(pop), "numeric")
})

test_that("discretePopSim clean", {
  pops<- list()
  pops$popMix<- discretePopSim(broods=1, b=1, j=.5, a=.5, breedFail=0,
                       varJ=0, varBreedFail=0, seasonVar=1,
                       sexRatio=NA, matingSystem=NA, N0=2, replicates=15, tf=10)
  pops$popExt<- discretePopSim(broods=1, b=1, j=.1, a=.2, breedFail=0,
                          varJ=0, varBreedFail=0, seasonVar=1,
                          sexRatio=NA, matingSystem=NA, N0=2, replicates=15, tf=10)
  pops$popMaxN<- discretePopSim(broods=1, b=10, j=.7, a=.9, breedFail=0,
                          varJ=0, varBreedFail=0, seasonVar=1,
                          sexRatio=NA, matingSystem=NA, N0=2, replicates=15, tf=10, maxN=1000)
  stats<- sapply(pops, summary)

  expect_equal(stats["decrease", "popExt"], 1)
  expect_equal(stats["increase", "popMaxN"], 1)
})
