context("Model")

test_that("combine LH and Env works", {
  expect_is(combineLH_Env(lh=LH(), env=Env()), "list")
  
  obj<- combineLH_Env(lh=LH(), env=Env())
  
  expect_equal(length(obj$seasonBroodEnv$broodEnv), nrow(obj$scenario))
  expect_equal(nrow(obj$seasonBroodEnv$parsBroodEnv), nrow(obj$scenario))
})

test_that("constructor works", {
  expect_is(Model(), "Model")
  expect_is(Model(lh=LH(), env=Env(), sim=Sim.discretePopSim()), "Model")
  expect_is(Model(sim=Sim.numericDistri()), "Model")
  expect_is(Model(sim=Sim.ssa()), "Model")

  ## A data.frame can't define a model. A Sim object is necessary.
  # obj<- Model()
  # expect_is(Model(S3Part(obj)), "Model")
  # expect_equivalent(Model(S3Part(obj)), obj)
})

test_that("subsetting works", {
  obj<- Model()
  expect_is(obj[c(1,4,8),], "Model")
  expect_is(head(obj), "Model")
  expect_is(obj[1:10,], "Model")

  expect_is(obj$lambda[2], "numeric")
  expect_is(obj[[2]], "numeric")
})
