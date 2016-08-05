context("Run models")

test_that("discrete time models works", {
  lh<- LH()
  sim<- Sim.discretePopSim()
  env<- Env()
  env<- Env(env[c(1,3,4,6),]) # stable environment only
  model<- Model(lh=lh, env=env, sim=sim)
  
  res<- run(model)
  resD<- result(res)
  plot(resD)
  resD<- result(res, type="Ntf")
})


test_that("compound distribution works", {
  lh<- LH()
  sim<- Sim.numericDistri()
  env<- Env()
  env<- Env(env[env$var == 0 & env$seasonAmplitude == 0,]) #stable environment
  model<- Model(lh=lh, env=env, sim=sim)
  
  res<- run(model)
  resD<- result(res)
  plot(resD)
  # TODO: fix wrong distributions!
  lapply(res@sim@raw, function(x) sum(x[[1]]$p))
  lapply(res@sim@raw, function(x) sum(x[[2]]$p))
})

test_that("IBM LH-behavior works", {
  tf <- 5 # Final time
  x0<- c(N1s=0, N1b=1, N1bF=0, N1j=0, N2s=0, N2b=1, N2bF=0, N2j=0)
  x0L<- lapply(c(2,10), function(x) x0 * x)
  
  params<- getParams.LH_Beh()
  
  transitionMat=transitionMat.LH_Beh; rateFunc=rateFunc.LH_Beh
  replicates<- 100
  cores<- 4
  
  sim<- Sim.ssa(N0=x0L, transitionMat=transitionMat.LH_Beh, rateFunc=rateFunc.LH_Beh, 
                tf=tf, replicates=replicates, raw=FALSE, Ntf=TRUE, stats=TRUE)
  model<- Model.ssa(pars=params, sim=sim)
  system.time(res<- run(model, cores=cores))
  
  resD<- result(res)
  plot(resD)
  resD<- result(res, type="Ntf")
})
