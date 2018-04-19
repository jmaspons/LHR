context("Class Simulation")


test_that("constructor", {
  expect_is(Sim(), "Sim")
  expect_is(Sim.discretePopSim(), "Sim")
  expect_is(Sim.discretePopSim(), "Sim.discretePopSim")
  expect_is(Sim.numericDistri(), "Sim")
  expect_is(Sim.numericDistri(), "Sim.numericDistri")
  expect_is(Sim.ABM(), "Sim")
  expect_is(Sim.ABM(), "Sim.ABM")

  expect_is(Sim(Model()), "Sim")
})

# test_that("subsetting", {
#   obj<- Sim()
#   expect_is(obj[c(1,4,8),], "Sim")
#   expect_is(head(obj), "Sim")
#   expect_is(obj[1:10,], "Sim")
#   
#   expect_is(obj$lambda[2], "numeric")
#   expect_is(obj[[2]], "numeric")
# })
