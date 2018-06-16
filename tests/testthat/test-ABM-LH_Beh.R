context("discreteABMSim LH_behavior")

test_that("discreteABMSim", {
  obj<- discreteABMSim()
  expect_is(obj, "discreteABMSim")
  obj<- discreteABMSim(Ntf=TRUE)
  expect_is(obj, "discreteABMSim")
  expect_is(summary(obj), "numeric")
  
  pop<- discreteABMSim2discretePopSim(obj)
  expect_is(pop, "discretePopSim")
  expect_is(summary(pop, dt = 1), "numeric")
  expect_is(r(pop, dt = 1), "matrix")
  expect_is(lambda(pop, dt = 1), "matrix")
  expect_is(trendsProp(pop, dt = 1), "numeric")
  
  expect_is(G(pop), "numeric")
  expect_is(Gmean(pop), "numeric")
})

test_that("sample parameter space", {
  
  expect_is(getScenario("identicalHab"), "numeric")
  expect_is(getScenario("mortalHab2"), "numeric")
  expect_is(getScenario("nestPredHab2"), "numeric")
  
  expect_is(setScenario(habDiffScenario="identicalHab"), "data.frame")
  expect_is(setScenario(habDiffScenario="mortalHab2"), "data.frame")
  expect_is(setScenario(habDiffScenario="nestPredHab2"), "data.frame")
  
  obj<- setParams2diff1(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5))
  expect_is(obj, "data.frame")
  
  expect_is(getParams.LH_Beh(), "data.frame")
  obj<- getParamsCombination.LH_Beh()
  expect_is(obj, "data.frame")
})

