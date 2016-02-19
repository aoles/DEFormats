library("DESeq2")
source("utils.R")

context("Simulate normalization factors")

test_that("output is a matrix of desired dimensions", {
  expect_is(simulateNormFactors(), "matrix")
  n = 100L
  expect_equal(nrow(simulateNormFactors(n = n)), n)
  m = 4L
  expect_equal(ncol(simulateNormFactors(m = m)), m)
})

test_that("the counts table is random and reproducible", {
  expect_identical(simulateNormFactors(), simulateNormFactors())
  expect_different(simulateNormFactors(seed = 0L), simulateNormFactors(seed = 1L))
})

test_that("calling the function doesn't change the global state of the random number generator", {
  expect_rng_unchanged(simulateNormFactors)
})
