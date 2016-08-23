context("discreteABMSim LH_behavior")

test_that("discreteABMSim", {
  obj<- discreteABMSim()
  
  
  pop<- discreteABMSim2discretePopSim(obj)
  expect_is(pop, "discretePopSim")
  expect_is(summary(pop, dt = 1), "data.frame")
  expect_is(r(pop, dt = 1), "matrix")
  expect_is(lambda(pop, dt = 1), "matrix")
  expect_is(trendsProp(pop, dt = 1), "data.frame")
  
  expect_is(G(pop), "numeric")
  expect_is(Gmean(pop), "numeric")
})

test_that("sample parameter space", {
  
  getScenario("identicalHab")
  getScenario("mortalHab2")
  getScenario("nestPredHab2")
  
  setScenario(scenario="identicalHab")
  setScenario(scenario="mortalHab2")
  setScenario(scenario="nestPredHab2")
  
  setParams2diff1(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,  a1=.1,ab1=.25,j1=.25,  a2=.1,ab2=.25,j2=.25, AFR=1, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1j=.5))
  getParams.LH_Beh()
  obj<- getParamsCombination.LH_Beh()
})

test_that("exploreABM", {
  matrix(rep(5, 8), nrow=4, ncol=8, dimnames=list(replicates=NULL, state=c("N1s", "N1b", "N1bF", "N1j", "N2s", "N2b", "N2bF", "N2j")))
  obj<- exploreABM()
})
