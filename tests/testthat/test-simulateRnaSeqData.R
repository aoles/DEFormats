library("DESeq2")
source("utils.R")

context("Simulate RNA-seq data")

se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
dds = DESeqDataSet(se, design = ~ condition)

test_that("output is of desired type", {
  expect_is(simulateRnaSeqData(output = "matrix"), "matrix")
  expect_is(simulateRnaSeqData(output = "RangedSummarizedExperiment"), "RangedSummarizedExperiment")
})

test_that("the number of replicates is correct", {
  rep = 5L
  expect_equal(ncol(simulateRnaSeqData(rep = rep)), 2*rep)
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
