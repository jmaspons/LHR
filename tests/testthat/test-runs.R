context("Run models")

# lh<- LH()
# env<- Env()
# env<- Env(env[env$var == 0 & env$seasonAmplitude == 0,]) #stable environment

test_that("discrete time models works", {
  sim<- Sim.discretePopSim()
  lh<- LH()
  lh<- lh[lh$lambda == 1,]
  model<- Model(lh=lh, env=Env(), sim=sim)
  
  res<- run(model, cores=parallel::detectCores())
  expect_is(res, "Model")
  tmp<- lapply(unlist(res@sim@raw, recursive=FALSE), expect_is, class="discretePopSim")
  
  expect_is(result(res), "data.frame")
  expect_is(result(res, type="Ntf"), "data.frame")
})


test_that("compound distribution works", {
  sim<- Sim.numericDistri()
  env<- Env()
  env<- env[env$var == 0,] ## Errors if var != 1
  lh<- LH()
  lh<- lh[lh$lambda == 1,]
  
  model<- Model(lh=lh, env=env, sim=sim)
  
  res<- run(model, cores=parallel::detectCores())
  expect_is(res, "Model")
  tmp<- lapply(unlist(res@sim@raw, recursive=FALSE), expect_is, class="numericDistri")
  
  expect_is(result(res), "data.frame")
  # Not available for numericDistri expect_is(result(res, type="Ntf"), "data.frame")
  
  # TODO: fix wrong distributions!
  sapply(res@sim@raw, function(x) sum(x[[1]]$p)) - 1
  sapply(res@sim@raw, function(x) sum(x[[2]]$p)) - 1
})


test_that("IBM LH-behavior works", {
  tf <- 5 # Final time
  x0<- c(N1s=0, N1b=1, N1bF=0, N1j=0, N2s=0, N2b=1, N2bF=0, N2j=0)
  x0L<- lapply(c(2,10), function(x) x0 * x)
  
  params<- LHR:::getParams.LH_Beh()
  
  transitionMat=LHR:::transitionMat.LH_Beh; rateFunc=LHR:::rateFunc.LH_Beh
  replicates<- 100
  cores<- parallel::detectCores()
  
  sim<- Sim.ssa(N0=x0L, transitionMat=transitionMat, rateFunc=rateFunc, 
                tf=tf, replicates=replicates, raw=FALSE, Ntf=TRUE, stats=TRUE)
  model<- Model.ssa(pars=params, sim=sim)
  system.time(res<- run(model, cores=cores, dt=0.5))
  
  expect_is(resD<- result(res), "data.frame")
  # plot(resD[,-1])
  expect_is(result(res, type="Ntf"), "data.frame")
})
