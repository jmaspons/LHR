context("discreteABMSim LH_behavior")

test_that("discreteABMSim", {
  expect_is(discreteABMSim(), "discreteABMSim")
  expect_is(popABM<- discreteABMSim(Ntf=TRUE), "discreteABMSim")
  expect_is(summary(popABM), "numeric")
  
  ## Test plots
  expect_is(plot(popABM), "NULL")
  expect_is(plotLH_behavior(popABM), "ggplot")
  
  pop<- discreteABMSim2discretePopSim(popABM)
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
  
  obj<- setParams2diff1(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,
                                          a1=.25,ab1=.2,sa1=.2,j1=.1,  a2=.25,ab2=.2,sa2=.2,j2=.1,
                                          AFR=2, K=500, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1sa=.5, P1j=.5))
  expect_is(obj, "data.frame")
  
  expect_is(getParams.LH_Beh(), "data.frame")
  expect_is(getParamsCombination.LH_Beh(), "data.frame")
})

test_that("plots", {
  model<- Model(pars=getParamsCombination.LH_Beh(habDiffScenario="nestPredHab2", behavior="skip"), type="ABM")
  res<- run(model)
  N0_Pest<- findN0_Pest(model)
  
  ## Test specific LH_behavior plots
  expect_is(plotLH_behavior(res, resultType="Pest_N0"), "ggplot")
  expect_is(plotLH_behavior(res, resultType="G"), "ggplot")
  expect_is(plotLH_behavior(res, resultType="Ntf"), "ggplot")
  
  expect_is(plotLH_behavior(N0_Pest, resultType="N0_Pest"), "ggplot")
  expect_is(plotLH_behavior(res@sim@raw[[1]][[1]], groups="all"), "ggplot")
  expect_is(plotLH_behavior(res@sim@raw[[1]][[1]], groups="habitat"), "ggplot")
  expect_is(plotLH_behavior(res@sim@raw[[1]][[1]], groups="age"), "ggplot")
  expect_is(plotLH_behavior(res@sim@raw[[1]][[1]], groups="habitat*age"), "ggplot")
  # plotLH_behavior(res@sim@raw$`LHslow-L1.05_Env1_nestPredHab2_skip`$`64`)
})
