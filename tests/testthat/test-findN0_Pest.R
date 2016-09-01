context("findN0_Pest")

test_that("discretePopSim", {
if (skip_on_cran()){
  sim<- Sim.discretePopSim(replicates=1000)
  lh<- LH(lambda=1.1, broods=2)[1:3,]
  env<- Env(seasonAmplitude=0, var=0, breedFail=.3)
  model<- Model(lh=lh, env=env, sim=sim)
  
  expect_is(findN0_Pest.scenario(scenario=data.frame(model)[1,], sim=sim, Pobjective=.5), "data.frame")
  expect_is(N0_Pest<- findN0_Pest(model=model, Pobjective=.5), "Model")
  expect_is(N0_Pest@sim@N0_Pest, "data.frame")
  expect_is(result(N0_Pest, type="N0_Pest"), "data.frame")
}
})

test_that("numericDistri", {
if (skip_on_cran()){
  sim<- Sim.numericDistri()
  lh<- LH(lambda=1.1, broods=2)[1:3,]
  env<- Env(seasonAmplitude=0, var=0, breedFail=.3)
  model<- Model(lh=lh, env=env, sim=sim)
  
  expect_is(findN0_Pest.scenario(scenario=data.frame(model)[1,], sim=sim, Pobjective=.5), "data.frame")
  expect_is(N0_Pest<- findN0_Pest(model=model, Pobjective=.5), "Model")
  expect_is(N0_Pest@sim@N0_Pest, "data.frame")
  expect_is(result(N0_Pest, type="N0_Pest"), "data.frame")
}
})

test_that("discreteABM", {
if (skip_on_cran()){
  sim<- Sim.ABM(replicates=1000)
  lh<- LH(lambda=1.1, broods=2)[1:3,]
  env<- Env(seasonAmplitude=0, var=0, breedFail=.3)
  model<- Model(lh=lh, env=env, sim=sim, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
  
  expect_is(findN0_Pest.scenario(scenario=data.frame(model)[1,], sim=sim, Pobjective=.5), "data.frame")
  expect_is(N0_Pest<- findN0_Pest(model=model, Pobjective=.5), "Model")
  expect_is(N0_Pest@sim@N0_Pest, "data.frame")
  expect_is(result(N0_Pest, type="N0_Pest"), "data.frame")
}
})


## Not ready
# test_that("ssa", {
# if (skip_on_cran()){
#   sim<- Sim.ssa(replicates=1000)
#   lh<- LH(lambda=1.1, broods=2)[1,]
#   env<- Env(seasonAmplitude=0, var=0, breedFail=.3)
#   model<- Model(lh=lh, env=env, sim=sim, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
#   
#   expect_is(findN0_Pest.scenario(scenario=data.frame(model)[1,], sim=sim, Pobjective=.5), "data.frame")
#   expect_is(N0_Pest<- findN0_Pest(model=model, Pobjective=.5), "Model")
#   expect_is(N0_Pest@sim@N0_Pest, "data.frame")
#   expect_is(result(N0_Pest, type="N0_Pest"), "data.frame")
# }
# })

