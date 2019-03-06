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

test_that("discreteABMSim extreme values", {
  ## Extinction
  suppressWarnings( lh<- LH(lambda=.3, broods=1, a=.3, AFR=c(1, 3), method="regular") )
  env<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.3)
  pars<- getParamsCombination.LH_Beh(lh, env, habDiffScenario="mortalHab2", behavior="preferHab2")
  
  N0<- c(N1s=0, N1b=1, N1bF=0, N2s=0, N2b=1, N2bF=0)
  sim<- Sim.ABM(N0=N0, replicates=1000, maxN=1000, raw=TRUE, discretePopSim=TRUE)
  model<- Model(pars=pars, sim=sim)
  
  resExt<- run(model)
  
  stats<- S3Part(resExt@sim)
  Ntf<- resExt@sim@Ntf
  pop<- resExt@sim@discretePopSim
  popABM<- resExt@sim@raw
  
  expect_equal(unique(stats$extinct), 1)
  expect_equal(unique(stats$decrease), 1)
  expect_equal(unique(stats$stable), 0)
  expect_equal(unique(stats$increase), 0)
  
  expect_equal(sum(Ntf[, ncol(Ntf)]), 0) # Last column in Ntf is the replicate with larger N.
  tmp<- lapply(unlist(pop, recursive=FALSE), function(x) expect_equal(unique(x[, ncol(x)]), NA_real_))
  tmp<- lapply(unlist(popABM, recursive=FALSE), function(x) expect_equal(sum(unique(x[,, dim(x)[3]])), 0))
  
  # maxN
  lh<- LH(lambda=2, broods=2, a=.8, AFR=c(1, 3), method="regular")
  env<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=.3)
  pars<- getParamsCombination.LH_Beh(lh, env, habDiffScenario="identicalHab", behavior="neutral")
  
  N0<- c(N1s=0, N1b=100, N1bF=0, N2s=0, N2b=100, N2bF=0)
  sim<- Sim.ABM(N0=N0, replicates=1000, maxN=1000, tf=100, raw=TRUE, discretePopSim=TRUE)
  model<- Model(pars=pars, sim=sim)
  
  resMaxN<- run(model)
  
  stats<- S3Part(resMaxN@sim)
  Ntf<- resMaxN@sim@Ntf
  pop<- resMaxN@sim@discretePopSim
  popABM<- resMaxN@sim@raw
  
  expect_equal(unique(stats$extinct), 0)
  expect_equal(unique(stats$decrease), 0)
  expect_equal(unique(stats$stable), 0)
  expect_equal(unique(stats$increase), 1)
  
  ## check model-discreteABMSim.R stop conditions (break).  It fills all classes (dim(popABM[[1]][[1]])[2]) with maxN
  ## Model.R -> runScenario.ABM() -> if (pars$Ntf) take the last N before all classes are MaxN
  expect_gt(sum(Ntf[, ncol(Ntf)]), sum(nrow(Ntf) * resMaxN@sim@params$maxN))
  tmp<- lapply(unlist(pop, recursive=FALSE), function(x) expect_equal(unique(x[, ncol(x)]), NA_real_))
  tmp<- lapply(unlist(popABM, recursive=FALSE), function(x) expect_true(all(is.na(unique(x[,, dim(x)[3]])))))
  
  # sapply(unlist(popABM, recursive=FALSE), summary)
})

test_that("sample parameter space", {
  expect_is(getDiffHabScenario("identicalHab"), "numeric")
  expect_is(getDiffHabScenario("mortalHab2"), "numeric")
  expect_is(getDiffHabScenario("nestPredHab2"), "numeric")
  
  expect_is(setHabScenario(habDiffScenario="identicalHab"), "data.frame")
  expect_is(setHabScenario(habDiffScenario="mortalHab2"), "data.frame")
  expect_is(setHabScenario(habDiffScenario="nestPredHab2"), "data.frame")
  
  obj<- getParams2diff1(params=data.frame(b1=1, b2=1,   broods=1, PbF1=.4, PbF2=.4,
                                          a1=.25,ab1=.2,sa1=.2,j1=.1,  a2=.25,ab2=.2,sa2=.2,j2=.1,
                                          AFR=2, Pb1=1, Pb2=1, c1=1, c2=1, cF=1, P1s=.5, P1b=.5, P1sa=.5, P1j=.5))
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
