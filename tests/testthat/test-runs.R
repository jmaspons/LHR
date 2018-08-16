context("Run discretePopSim models")

test_that("discrete time models", {
  lh<- LH(method="LH axes")
  env<- Env()
  
  ## Female only
  sim<- Sim.discretePopSim()
  model<- Model(lh=lh, env=env, sim=sim)
  
  if (skip_on_cran()){
    res<- run(model)
    expect_is(res, "Model")
    tmp<- lapply(unlist(res@sim@raw, recursive=FALSE), expect_is, class=c("discretePopSim", "logical"))
  
    expect_is(result(res), "data.frame")
    expect_is(result(res, type="Ntf"), "data.frame")
    
    ## Test subsetting
    expect_identical(length(res@sim@raw), nrow(res))
    expect_identical(length(res[1,]@sim@raw), nrow(res[1,]))
    expect_identical(length(res[c(1,3),]@sim@raw), nrow(res[c(1,3),]))
    
    expect_equal(nrow(res@sim@Ntf) / length(res@sim@params$N0), nrow(res))
    expect_equal(nrow(res[1,]@sim@Ntf) / length(res@sim@params$N0), nrow(res[1,]))
    expect_equal(nrow(res[c(1,3),]@sim@Ntf) / length(res@sim@params$N0), nrow(res[c(1,3),]))
    
    ## Test rbind
    expect_is(rbind(res[1:3,], res[4:6,]), "Model")
    expect_identical(rbind(res[1:3,], res[4:6,]), res[1:6,])
    ## Duplicated ids
    lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
    lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
    env1<- Env(seasonAmplitude=0, varJ=0, varA=0)
    env2<- Env(varJ=0, varA=0, breedFail=0.5)
    model1<- Model(lh=lh1, env=env1, sim=sim)
    model2<- Model(lh=lh2, env=env2, sim=sim)
    res1<- run(model1)
    res2<- run(model2)
    res12<- rbind(res1, res2)
    expect_identical(nrow(res12), nrow(res1) + nrow(res2))
    expect_setequal(res12@sim@Ntf$idScenario, res12$idScenario)
    expect_setequal(names(res12@sim@discretePopSim), res12$idScenario)
    
    ## Test plots
    expect_is(plot(res, resultType="Pest_N0"), "ggplot")
    expect_is(plot(res, resultType="G"), "ggplot")
    expect_equal(plot(res, resultType="N0_Pest"), NA)
    expect_is(plot(res, resultType="Ntf"), "ggplot")
    
    expect_is(hist(res, resultType="Ntf"), "ggplot")
  }
})

test_that("discrete models with 2 sexes", {
  lh<- LH(method="LH axes")
  env<- Env()
  
  ## 2 sexes
  sim<- Sim.discretePopSim(sexRatio=0.5, matingSystem="monogamy")
  model<- Model(lh=lh, env=env, sim=sim)
  
  if (skip_on_cran()){
    res<- run(model)
    expect_is(res, "Model")
    
    popList<- unlist(res@sim@raw, recursive=FALSE)
    
    # TODO: Not implemented models return NA
    popList<- popList[sapply(popList, function(x) !all(is.na(x)))]
    
    tmp<- lapply(popList, expect_is, class="discretePopSim")
    
    ## Test subsetting
    expect_identical(length(res@sim@raw), nrow(res))
    expect_identical(length(res[1,]@sim@raw), nrow(res[1,]))
    expect_identical(length(res[c(1,3),]@sim@raw), nrow(res[c(1,3),]))
    
    expect_equal(nrow(res@sim@Ntf) / length(res@sim@params$N0), nrow(res))
    expect_equal(nrow(res[1,]@sim@Ntf) / length(res@sim@params$N0), nrow(res[1,]))
    expect_equal(nrow(res[c(1,3),]@sim@Ntf) / length(res@sim@params$N0), nrow(res[c(1,3),]))
    
    ## Test rbind
    expect_is(rbind(res[1:3,], res[4:6,]), "Model")
    expect_identical(rbind(res[1:3,], res[4:6,]), res[1:6,])
    expect_error(rbind(res, res)) # duplicated scenarios
    ## Duplicated ids
    lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
    lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
    env1<- Env(seasonAmplitude=0, varJ=0, varA=0)
    env2<- Env(varJ=0, varA=0, breedFail=0.5)
    model1<- Model(lh=lh1, env=env1, sim=sim)
    model2<- Model(lh=lh2, env=env2, sim=sim)
    res1<- run(model1)
    res2<- run(model2)
    res12<- rbind(res1, res2)
    expect_identical(nrow(res12), nrow(res1) + nrow(res2))
    expect_setequal(res12@sim@Ntf$idScenario, res12$idScenario)
    expect_setequal(names(res12@sim@discretePopSim), res12$idScenario)
    
    ## Test plots
    expect_is(plot(res, resultType="Pest_N0"), "ggplot")
    expect_is(plot(res, resultType="G"), "ggplot")
    expect_equal(plot(res, resultType="N0_Pest"), NA)
    expect_is(plot(res, resultType="Ntf"), "ggplot")
    
    expect_is(hist(res, resultType="Ntf"), "ggplot")
  }
})


context("Run numericDistri models")

test_that("compound distribution", {
  sim<- Sim.numericDistri()
  env<- Env(varJ=0, varA=0) ## Errors if var != 0
  lh<- LH(method="LH axes")
  
  model<- Model(lh=lh, env=env, sim=sim)
  
  if (skip_on_cran()){
    res<- run(model)
    expect_is(res, "Model")
    
    distriList<- unlist(res@sim@raw, recursive=FALSE)
    
    ## TODO: Not implemented models return NA
    distriList<- distriList[sapply(distriList, function(x) !all(is.na(x)))]
    
    tmp<- lapply(distriList, expect_is, class="numericDistri")
    
    expect_is(result(res), "data.frame")
    # Not available for numericDistri expect_is(result(res, type="Ntf"), "data.frame")
    
    tmp<- lapply(distriList, function(x) expect_gt(abs(sum(x$p)), 0.95))
    
    ## Test subsetting
    expect_identical(length(res@sim@raw), nrow(res))
    expect_identical(length(res[1,]@sim@raw), nrow(res[1,]))
    expect_identical(length(res[c(1,3),]@sim@raw), nrow(res[c(1,3),]))
    
    # Ntf<- result(res, type="Ntf")
    expect_equal(nrow(result(res, type="Ntf")) / length(res@sim@params$N0), nrow(res))
    expect_equal(nrow(result(res[1,], type="Ntf")) / length(res@sim@params$N0), nrow(res[1,]))
    expect_equal(nrow(result(res[c(1,3),], type="Ntf")) / length(res@sim@params$N0), nrow(res[c(1,3),]))
    
    ## Test rbind
    expect_is(rbind(res[1:3,], res[4:6,]), "Model")
    expect_identical(rbind(res[1:3,], res[4:6,]), res[1:6,])
    expect_error(rbind(res, res)) # duplicated scenarios
    ## Duplicated ids
    lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
    lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
    env1<- Env(seasonAmplitude=0, varJ=0, varA=0)
    env2<- Env(varJ=0, varA=0, breedFail=0.5)
    model1<- Model(lh=lh1, env=env1, sim=sim)
    model2<- Model(lh=lh2, env=env2, sim=sim)
    res1<- run(model1)
    res2<- run(model2)
    res12<- rbind(res1, res2)
    expect_identical(nrow(res12), nrow(res1) + nrow(res2))
    expect_setequal(names(res12@sim@raw), res12$idScenario)

    ## Test plots
    expect_is(plot(res, resultType="Pest_N0"), "ggplot")
    expect_is(plot(res, resultType="G"), "ggplot")
    expect_equal(plot(res, resultType="N0_Pest"), NA)
    expect_is(plot(res, resultType="Ntf"), "ggplot")
  }
})

test_that("compound distribution with environmental variation", {
  sim<- Sim.numericDistri()
  env<- Env()
  env<- env[env$varJ != 0,] ## Errors if var != 0
  lh<- LH(method="LH axes")
  
  model<- Model(lh=lh, env=env, sim=sim)
  
  if (skip_on_cran()){
    res<- run(model)
    expect_is(res, "Model")
    
    distriList<- unlist(res@sim@raw, recursive=FALSE)
    
    # TODO: Not implemented models return NA
    distriList<- distriList[sapply(distriList, function(x) !all(is.na(x)))]
    
    tmp<- lapply(distriList, expect_is, class="numericDistri")
    
    expect_is(result(res), "data.frame")
    # Not available for numericDistri expect_is(result(res, type="Ntf"), "data.frame")
    
    # TODO: fix wrong distributions!
    ## Fails for seasonAmplitude=1 & var=0.1 and when breedFail=1
    distriList<- distriList[sapply(distriList, function(x) !all(is.na(x$p)))] ## TODO: fix it! remove some results where probability is NA
    tmp<- sapply(distriList, function(x) expect_gt(abs(sum(x$p)), 0.95))
    
    ## Test subsetting
    expect_identical(length(res@sim@raw), nrow(res))
    expect_identical(length(res[1,]@sim@raw), nrow(res[1,]))
    expect_identical(length(res[c(1,3),]@sim@raw), nrow(res[c(1,3),]))
    
    # Ntf<- result(res, type="Ntf")
    expect_equal(nrow(result(res, type="Ntf")) / length(res@sim@params$N0), nrow(res))
    expect_equal(nrow(result(res[1,], type="Ntf")) / length(res@sim@params$N0), nrow(res[1,]))
    expect_equal(nrow(result(res[c(1,3),], type="Ntf")) / length(res@sim@params$N0), nrow(res[c(1,3),]))
    
    ## Test rbind
    expect_is(rbind(res[1:3,], res[4:6,]), "Model")
    expect_identical(rbind(res[1:3,], res[4:6,]), res[1:6,])
    expect_error(rbind(res, res)) # duplicated scenarios
    ## Duplicated ids
    lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
    lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
    env1<- Env(seasonAmplitude=0, varJ=0, varA=.01)
    env2<- Env(varJ=0, varA=0.05, breedFail=0.5)
    model1<- Model(lh=lh1, env=env1, sim=sim)
    model2<- Model(lh=lh2, env=env2, sim=sim)
    res1<- run(model1)
    res2<- run(model2)
    res12<- rbind(res1, res2)
    expect_identical(nrow(res12), nrow(res1) + nrow(res2))
    expect_setequal(names(res12@sim@raw), res12$idScenario)
    
    ## Test plots
    expect_is(plot(res, resultType="Pest_N0"), "ggplot")
    expect_is(plot(res, resultType="G"), "ggplot")
    expect_equal(plot(res, resultType="N0_Pest"), NA)
    expect_is(plot(res, resultType="Ntf"), "ggplot")
  }
})

context("Run discreteABMSim models")

test_that("ABM LH-behavior", { 
  lh<- LH(method="LH axes")
  env<- Env(seasonAmplitude=0, varJ=0, varA=0)
  sim<- Sim.ABM()
  pars<- getParamsCombination.LH_Beh(lh=lh, env=env, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
  model<- Model(sim=sim, pars=pars)
  
  ## TODO: make it work
  # model<- Model(lh=lh, env=env, sim=sim)
  # model<- model[model$habDiff == "nestPredHab2" & model$behavior == "learnExploreBreed", ]
  
  if (skip_on_cran()){
    res<- run(model)
    expect_is(res, "Model")
    
    popABML<- unlist(res@sim@raw, recursive=FALSE)
    tmp<- lapply(popABML, expect_is, class="discreteABMSim")
    
    expect_is(result(res), "data.frame")
    expect_is(result(res, type="Ntf"), "data.frame")
    
    ## Test subsetting
    expect_identical(length(res@sim@raw), nrow(res))
    expect_identical(length(res[1,]@sim@raw), nrow(res[1,]))
    expect_identical(length(res[c(1,3),]@sim@raw), nrow(res[c(1,3),]))
    
    expect_equal(nrow(res@sim@Ntf) / length(res@sim@params$N0), nrow(res))
    expect_equal(nrow(res[1,]@sim@Ntf) / length(res@sim@params$N0), nrow(res[1,]))
    expect_equal(nrow(res[c(1,3),]@sim@Ntf) / length(res@sim@params$N0), nrow(res[c(1,3),]))
    
    ## Test rbind
    expect_is(rbind(res[1:3,], res[4:6,]), "Model")
    expect_identical(rbind(res[1:3,], res[4:6,]), res[1:6,])
    # TODO: differences in row sorting
    # a<- res[1:6,]
    # b<- rbind(res[1:3,], res[4:6,])
    # all.equal(a@sim@Ntf[order(a@sim@Ntf$idScenario),], b@sim@Ntf[order(b@sim@Ntf$idScenario),])
    expect_error(rbind(res, res)) # duplicated scenarios

    ## Duplicated ids
    lh1<- LH(lambda=1, broods=1, a=.7, method="regular")
    lh2<- LH(lambda=1.1, broods=1, a=.6, method="regular")
    env1<- Env(seasonAmplitude=0, varJ=0, varA=0)
    env2<- Env(seasonAmplitude=0, varJ=0, varA=0, breedFail=0.4)
    model1<- Model(lh=lh1, env=env1, sim=sim, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
    model2<- Model(lh=lh2, env=env2, sim=sim, habDiffScenario="nestPredHab2", behavior="learnExploreBreed")
    res1<- run(model1)
    res2<- run(model2)
    res12<- rbind(res1, res2)
    expect_identical(nrow(res12), nrow(res1) + nrow(res2))
    expect_setequal(res12@sim@Ntf$idScenario, res12$idScenario)
    expect_setequal(names(res12@sim@discretePopSim), res12$idScenario)
    expect_setequal(names(res12@sim@raw), res12$idScenario)
    
    ## Test plots
    expect_is(plot(res, resultType="Pest_N0"), "ggplot")
    expect_is(plot(res, resultType="G"), "ggplot")
    expect_equal(plot(res, resultType="N0_Pest"), NA)
    expect_is(plot(res, resultType="Ntf"), "ggplot")
    
    expect_is(hist(res, resultType="Ntf"), "ggplot")
  }
})


