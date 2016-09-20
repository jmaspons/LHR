context("Class Environment")

test_that("constructor", {
  expect_is(Env(), "Env")
  expect_is(Env(varJ=.1), "Env")
  expect_is(Env(seasonAmplitude=.3), "Env")
  expect_is(Env(seasonRange=matrix(c(0,1, .4,1, .8,1, 1, 1), 2, 2, byrow=TRUE)), "Env") ## ERROR
  expect_is(Env(varJ=.1, seasonAmplitude=.3), "Env")
  expect_is(Env(varJ=1), "Env") # parameters out of the Beta distribution domain
  expect_error(Env(seasonRange=c(0,1)), "Env") # mean and range are redundant parameters on a sinusoidal function
  
  obj<- Env()
  expect_is(Env(S3Part(obj)), "Env")
  expect_equivalent(Env(S3Part(obj)), obj)
})

test_that("subsetting", {
  obj<- Env()
  expect_is(obj[c(1,4,8),], "Env")
  expect_is(head(obj), "Env")
  expect_is(obj[1:10,], "Env")
  
  expect_is(obj$seasonAmplitude, "numeric")
  expect_is(obj[[2]], "numeric")
})

test_that("seasonal pattern", {
  env<- Env()
  expect_is(seasonalPattern(env), "matrix")
  expect_is(seasonalPattern(env, resolution=365), "matrix")
  expect_is(seasonalPattern(env, cicles=2), "matrix")
  
  expect_is(seasonOptimCal(env), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=1, criterion="maxFirst"), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=3, criterion="maxFirst"), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=1, criterion="maxMean"), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=3, criterion="maxMean"), "list")
  
  
  env<- data.frame(Env())
  expect_is(seasonalPattern(env), "matrix")
  expect_is(seasonalPattern(env, resolution=365), "matrix")
  expect_is(seasonalPattern(env, cicles=2), "matrix")
  
  expect_is(seasonOptimCal(env), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=1, criterion="maxFirst"), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=3, criterion="maxFirst"), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=1, criterion="maxMean"), "list")
  expect_is(seasonOptimCal(env, resolution=12, nSteps=3, interval=3, criterion="maxMean"), "list")
})
