context("Beta (Negative) Binomial distribution")

test_that("betabinom works", {
  expect_is(dbetabinom(2, 5, 1:9,9:1), "numeric")
  expect_is(pbetabinom(2, 5, 2,3), "numeric")
  expect_is(qbetabinom(.6, 5, 2,3), "numeric")
  expect_is(rbetabinom(10, 5, 2,3), "numeric")
})

test_that("betanbinom works", {
  expect_is(dbetanbinom(2, 5, 1:9,9:1), "numeric")
  expect_is(pbetanbinom(2, 5, 2,3), "numeric")
  expect_is(qbetanbinom(.6, 5, 2,3), "numeric")
  expect_is(rbetanbinom(10, 5, 2,3), "numeric")
})
