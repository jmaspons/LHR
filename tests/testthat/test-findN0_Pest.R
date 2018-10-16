context("findN0_Pest")

## discretePopSim ----
test_that("discretePopSim", {
if (skip_on_cran()){
  sim<- Sim.discretePopSim(replicates=1000)
  lh<- LH(lambda=1.1, broods=2)[1:3,]
  env<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.3)
  model<- Model(lh=lh, env=env, sim=sim)
  
  expect_is(findN0_Pest.scenario(scenario=data.frame(model)[1,], sim=sim, Pobjective=.5), "data.frame")
  expect_is(N0_Pest<- findN0_Pest(model=model, Pobjective=.5), "Model")
  expect_is(N0_Pest@sim@N0_Pest, "data.frame")
  expect_is(result(N0_Pest, type="N0_Pest"), "data.frame")
  
  ## Test subsetting
  expect_identical(nrow(N0_Pest@sim@N0_Pest), 3L)
  expect_identical(nrow(N0_Pest[1,]@sim@N0_Pest), 1L)
  expect_identical(nrow(N0_Pest[c(1,3),]@sim@N0_Pest), 2L)
  
  ## Test rbind
  lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
  lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
  env1<- Env(varJ=0, varA=0, breedFail=.3)
  env2<- Env(seasonAmplitude=0, varJ=0, varA=0)
  model1<- Model(lh=lh1, env=env1, sim=sim)
  model2<- Model(lh=lh2, env=env2, sim=sim)
  N0_Pest1<- findN0_Pest(model1)
  N0_Pest2<- findN0_Pest(model2)
  N0_Pest12<- rbind(N0_Pest1, N0_Pest2)
  expect_identical(nrow(N0_Pest12), nrow(N0_Pest1) + nrow(N0_Pest2))
  expect_setequal(N0_Pest12@sim@N0_Pest$idScenario, N0_Pest12$idScenario)
  expect_identical(rownames(N0_Pest12), N0_Pest12$idScenario)
  
  ## Test plots
  expect_equal(plot(N0_Pest, resultType="Pest_N0"), NA)
  expect_equal(plot(N0_Pest, resultType="G"), NA)
  expect_is(plot(N0_Pest, resultType="N0_Pest"), "ggplot")
  expect_equal(plot(N0_Pest, resultType="Ntf"), NA)
  
  ## Critical values
  lh<- LH(lambda=.3, broods=1, a=.3, method="regular")
  # Pest < Pobjective for N0 == maxN
  sim<- Sim.discretePopSim(replicates=1000, maxN=1000)
  model<- Model(lh=lh, env=env, sim=sim)
  N0_Pest<- findN0_Pest(model=model, Pobjective=.5)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0_Pest), sim@params$maxN)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0interpoled), NA_real_)
  # Pest == 0 for N0 == maxN
  sim<- Sim.discretePopSim(replicates=1000, maxN=1000, tf=100)
  model<- Model(lh=lh, env=env, sim=sim)
  N0_Pest<- findN0_Pest(model=model, Pobjective=.5)
  expect_equal(unique(N0_Pest@sim@N0_Pest$Pest), 0)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0_Pest), sim@params$maxN)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0interpoled), NA_real_)
  # Pest > Pobjective for N0 == 1
  sim<- Sim.discretePopSim(replicates=1000, maxN=1000)
  model<- Model(lh=LH(), env=env, sim=sim)
  N0_Pest<- findN0_Pest(model=model, Pobjective=.1)
  tmp<- sapply(N0_Pest@sim@N0_Pest$Pest, expect_gt, 0.1)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0_Pest), 1)
  tmp<- sapply(N0_Pest@sim@N0_Pest$N0interpoled, expect_lt, 1)
}
})

## numericDistri ----
test_that("numericDistri", {
if (skip_on_cran()){
  sim<- Sim.numericDistri()
  lh<- LH(lambda=1.1, broods=2)[1:3,]
  env<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.3)
  model<- Model(lh=lh, env=env, sim=sim)
  
  expect_is(findN0_Pest.scenario(scenario=data.frame(model)[1,], sim=sim, Pobjective=.5), "data.frame")
  expect_is(N0_Pest<- findN0_Pest(model=model, Pobjective=.5), "Model")
  expect_is(N0_Pest@sim@N0_Pest, "data.frame")
  expect_is(result(N0_Pest, type="N0_Pest"), "data.frame")
  
  ## Test subsetting
  expect_identical(nrow(N0_Pest@sim@N0_Pest), 3L)
  expect_identical(nrow(N0_Pest[1,]@sim@N0_Pest), 1L)
  expect_identical(nrow(N0_Pest[c(1,3),]@sim@N0_Pest), 2L)
  
  ## Test rbind
  lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
  lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
  env1<- Env(varJ=0, varA=0, breedFail=.3)
  env2<- Env(seasonAmplitude=0, varJ=0, varA=0)
  model1<- Model(lh=lh1, env=env1, sim=sim)
  model2<- Model(lh=lh2, env=env2, sim=sim)
  N0_Pest1<- findN0_Pest(model1)
  N0_Pest2<- findN0_Pest(model2)
  N0_Pest12<- rbind(N0_Pest1, N0_Pest2)
  expect_identical(nrow(N0_Pest12), nrow(N0_Pest1) + nrow(N0_Pest2))
  expect_setequal(N0_Pest12@sim@N0_Pest$idScenario, N0_Pest12$idScenario)
  expect_identical(rownames(N0_Pest12), N0_Pest12$idScenario)
  
  ## Test plots
  expect_equal(plot(N0_Pest, resultType="Pest_N0"), NA)
  expect_equal(plot(N0_Pest, resultType="G"), NA)
  expect_is(plot(N0_Pest, resultType="N0_Pest"), "ggplot")
  expect_equal(plot(N0_Pest, resultType="Ntf"), NA)
  
  ## TODO: Critical values
}
})

## discreteABM ----
test_that("discreteABM", {
if (skip_on_cran()){
  sim<- Sim.ABM(replicates=1000)
  lh<- LH(lambda=1.1, broods=2)[1:3,]
  env<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.3)
  model<- Model(lh=lh, env=env, sim=sim, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
  
  expect_is(findN0_Pest.scenario(scenario=data.frame(model)[1,], sim=sim, Pobjective=.5), "data.frame")
  expect_is(N0_Pest<- findN0_Pest(model=model, Pobjective=.5), "Model")
  expect_is(N0_Pest@sim@N0_Pest, "data.frame")
  expect_is(result(N0_Pest, type="N0_Pest"), "data.frame")
  
  ## Test subsetting
  expect_identical(nrow(N0_Pest@sim@N0_Pest), 3L)
  expect_identical(nrow(N0_Pest[1,]@sim@N0_Pest), 1L)
  expect_identical(nrow(N0_Pest[c(1,3),]@sim@N0_Pest), 2L)
  
  ## Test rbind
  lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
  lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
  env1<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.3)
  env2<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.5)
  model1<- Model(lh=lh1, env=env1, sim=sim, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
  model2<- Model(lh=lh2, env=env2, sim=sim, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
  N0_Pest1<- findN0_Pest(model1)
  N0_Pest2<- findN0_Pest(model2)
  N0_Pest12<- rbind(N0_Pest1, N0_Pest2)
  expect_identical(nrow(N0_Pest12), nrow(N0_Pest1) + nrow(N0_Pest2))
  expect_setequal(N0_Pest12@sim@N0_Pest$idScenario, N0_Pest12$idScenario)
  expect_identical(rownames(N0_Pest12), N0_Pest12$idScenario)
  
  ## Test plots
  expect_equal(plot(N0_Pest, resultType="Pest_N0"), NA)
  expect_equal(plot(N0_Pest, resultType="G"), NA)
  expect_is(plot(N0_Pest, resultType="N0_Pest"), "ggplot")
  expect_equal(plot(N0_Pest, resultType="Ntf"), NA)
  
  ## Critical values
  lh<- LH(lambda=.3, broods=1, a=.3, method="regular")
  env<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.3)
  pars<- getParamsCombination.LH_Beh(lh, env, habDiffScenario="mortalHab2", behavior="preferHab2")
  # Pest < Pobjective for N0 == maxN
  sim<- Sim.ABM(replicates=1000, maxN=1000, raw=FALSE)
  model<- Model(pars=pars, sim=sim)
  N0_Pest<- findN0_Pest(model=model, Pobjective=.5)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0_Pest), sim@params$maxN)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0interpoled), NA_real_)
  # Pest == 0 for N0 == maxN
  sim<- Sim.ABM(replicates=1000, maxN=1000, tf=100, raw=FALSE)
  model<- Model(pars=pars, sim=sim)
  N0_Pest<- findN0_Pest(model=model, Pobjective=.5)
  expect_equal(unique(N0_Pest@sim@N0_Pest$Pest), 0)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0_Pest), sim@params$maxN)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0interpoled), NA_real_)
  # Pest > Pobjective for N0 == 1
  sim<- Sim.ABM(replicates=1000, maxN=1000, raw=FALSE)
  pars<- getParamsCombination.LH_Beh(LH(lambda=1.2), env, habDiffScenario="identicalHab", behavior="neutral")
  model<- Model(pars=pars, sim=sim)
  N0_Pest<- findN0_Pest(model=model, Pobjective=.1)
  tmp<- sapply(N0_Pest@sim@N0_Pest$Pest, expect_gt, 0.1)
  expect_equal(unique(N0_Pest@sim@N0_Pest$N0_Pest), 1)
  tmp<- sapply(N0_Pest@sim@N0_Pest$N0interpoled, expect_lt, 1)
}
})


