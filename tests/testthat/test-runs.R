context("Run discretePopSim models")

test_that("discrete time models", {
  lh<- LH()
  lh<- lh[lh$lambda == 1,]
  env<- Env()
  
  ## Female only
  sim<- Sim.discretePopSim()
  model<- Model(lh=lh, env=env, sim=sim)
  
  res<- run(model)
  expect_is(res, "Model")
  tmp<- lapply(unlist(res@sim@raw, recursive=FALSE), expect_is, class="discretePopSim")
  
  expect_is(result(res), "data.frame")
  expect_is(result(res, type="Ntf"), "data.frame")
})

test_that("discrete models with 2 sexes", {
  lh<- LH()
  lh<- lh[lh$lambda == 1,]
  env<- Env()
  
  ## 2 sexes
  sim<- Sim.discretePopSim(sexRatio=0.5, matingSystem="monogamy")
  model<- Model(lh=lh, env=env, sim=sim)
  
  res<- run(model)
  expect_is(res, "Model")
  
  popList<- unlist(res@sim@raw, recursive=FALSE)
  # Not implemented models return NA
  popList<- popList[sapply(popList, function(x) !all(is.na(x)))]
  
  tmp<- lapply(popList, expect_is, class="discretePopSim")
})


context("Run numericDistri models")

test_that("compound distribution", {
  sim<- Sim.numericDistri()
  env<- Env()
  env<- env[env$var == 0,] ## Errors if var != 0
  lh<- LH()
  lh<- lh[lh$lambda == 1,]
  
  model<- Model(lh=lh, env=env, sim=sim)
  
  res<- run(model)
  expect_is(res, "Model")
  
  distriList<- unlist(res@sim@raw, recursive=FALSE)
  # Not implemented models return NA
  distriList<- distriList[sapply(distriList, function(x) !all(is.na(x)))]
  tmp<- lapply(distriList, expect_is, class="numericDistri")
  
  expect_is(result(res), "data.frame")
  # Not available for numericDistri expect_is(result(res, type="Ntf"), "data.frame")
  
  # TODO: fix wrong distributions!
  tmp<- lapply(distriList, function(x) expect_gt(abs(sum(x$p)), 0.95))
})

test_that("compound distribution with environmental variation", {
  sim<- Sim.numericDistri()
  env<- Env()
  env<- env[env$var != 0,] ## Errors if var != 0
  lh<- LH()
  lh<- lh[lh$lambda == 1,]
  
  model<- Model(lh=lh, env=env, sim=sim)
  
  res<- run(model)
  expect_is(res, "Model")
  distriList<- unlist(res@sim@raw, recursive=FALSE)
  # Not implemented models return NA
  distriList<- distriList[sapply(distriList, function(x) !all(is.na(x)))]
  
  tmp<- lapply(distriList, expect_is, class="numericDistri")
  
  expect_is(result(res), "data.frame")
  # Not available for numericDistri expect_is(result(res, type="Ntf"), "data.frame")
  
  # TODO: fix wrong distributions!
  tmp<- sapply(distriList, function(x) expect_gt(abs(sum(x$p)), 0.95))
})

context("Run discreteABMSim models")

test_that("ABM LH-behavior", { 
  lh<- LH()
  lh<- lh[lh$lambda == 1,]
  env<- Env(seasonAmplitude=0, var=0)
  sim<- Sim.ABM()
  pars<- getParamsCombination.LH_Beh(lh=lh, env=env, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
  model<- Model(sim=sim, pars=pars)
  
  ## TODO: make it work
  # model<- Model(lh=lh, env=env, sim=sim)
  # model<- model[model$habDiff == "nestPredHab2" & model$behavior == "learnExploreBreed", ]
  
  res<- run(model)
  expect_is(res, "Model")
  
  #TODO:
  popABML<- res@sim@raw
  tmp<- lapply(popABML, expect_is, class="exploreABMSim")
  
  expect_is(result(res), "data.frame")
  expect_is(result(res, type="Ntf"), "data.frame")
})


context("Run ssa models")

test_that("IBM LH-behavior", {
  tf <- 2 # Final time
  replicates<- 10
  x0<- c(N1s=0, N1b=1, N1bF=0, N1j=0, N2s=0, N2b=1, N2bF=0, N2j=0)
  x0L<- lapply(c(2,10), function(x) x0 * x)
  
  params<- getParams.LH_Beh.ssa()
  transitionMat<- transitionMat.LH_Beh
  rateFunc<- rateFunc.LH_Beh
    
  sim<- Sim.ssa(N0=x0L, transitionMat=transitionMat, rateFunc=rateFunc, 
                tf=tf, replicates=replicates, raw=FALSE, Ntf=TRUE, stats=TRUE)
  model<- Model(pars=params, sim=sim)

  system.time(res<- run(model, dt=0.5))
  
  expect_is(resD<- result(res), "data.frame")
  # plot(resD[,-1])
  expect_is(result(res, type="Ntf"), "data.frame")
})
