context("Model")

test_that("constructor works", {
  expect_is(Model(), "Model")
  expect_is(Model(lh=LH(), env=Env()), "Model")
})
