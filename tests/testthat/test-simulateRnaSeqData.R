library("DESeq2")
source("utils.R")

context("Simulate RNA-seq data")

test_that("output is of desired type", {
  for(type in c("matrix", "RangedSummarizedExperiment"))
    expect_is(simulateRnaSeqData(output = type), type)
})

test_that("the count matrix has desired dimensions", {
  n = 100L
  expect_equal(nrow(simulateRnaSeqData(n = n)), n)
  m = 4L
  expect_equal(ncol(simulateRnaSeqData(m = m)), m)
})

test_that("the counts table is random and reproducible", {
  expect_identical(simulateRnaSeqData(), simulateRnaSeqData())
  expect_different(simulateRnaSeqData(seed = 0L), simulateRnaSeqData(seed = 1L))
})

test_that("calling the function doesn't change the global state of the random number generator", {
  expect_rng_unchanged(simulateRnaSeqData)
})
