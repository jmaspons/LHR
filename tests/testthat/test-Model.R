context("Class Model")

test_that("combine LH and Env", {
  expect_is(obj<- combineLH_Env(lh=LH(), env=Env()), "list")
  
  expect_equal(length(obj$seasonBroodEnv$broodEnv), nrow(obj$scenario))
  expect_equal(nrow(obj$seasonBroodEnv$parsBroodEnv), nrow(obj$scenario))
  
  d<- obj$scenario
  expect_true(all(abs(d$j - d$jbr * d$jind) < .Machine$double.eps))
})

test_that("Model class constructors", {
  lh<- LH()
  env<- Env()
  expect_is(Model(), "Model")
  expect_is(Model(lh=lh, env=Env(), sim=Sim.discretePopSim()), "Model")
  expect_is(Model(sim=Sim.numericDistri()), "Model")

  ## No seasonality implemented
  expect_is(Model(env=env[env$seasonAmplitude == 0,], sim=Sim.ABM(),
                  patchScenario=getPatchScenario(habDiffScenario="nestPredHab2", behavior="learnExploreBreed")), "Model")
})


test_that("Model subclasses constructors", {
  expect_is(Model(type="discretePopSim"), "Model.discretePopSim")
  expect_is(Model(type="numericDistri"), "Model.numericDistri")
  expect_is(Model(sim=Sim.ABM(), env=Env(seasonAmplitude=0)), "Model.ABM")
})

test_that("subsetting and rbind", {
  obj<- Model()
  expect_is(obj[c(1,4,8),], "Model")
  expect_is(head(obj), "Model")
  expect_is(obj[1:10,], "Model")

  expect_is(obj$lambda[2], "numeric")
  expect_is(obj[[1]], "character")
  
  expect_is(rbind(obj[1:3,], obj[4:6,]), "Model")
  expect_error(rbind(obj, obj)) # duplicated scenarios
  expect_identical(rbind(obj[1:3,], obj[4:6,]), obj[1:6,])
})

test_that("parameter space plots", {
  expect_is(plot(Model(type="discretePopSim")), "NULL")
  expect_is(plot(Model(type="numericDistri")), "NULL")
  expect_is(plot(Model(sim=Sim.ABM(), env=Env(seasonAmplitude=0))), "NULL")
})
