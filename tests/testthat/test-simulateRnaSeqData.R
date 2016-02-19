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

test_that("columns are named after conditions", {
  conditions = c("a", "b")
  n = colnames(simulateRnaSeqData(conditions = conditions))
  n = sub(" [0-9]+$", "",  n)
  expect_identical(n , rep(conditions, each = length(n)/length(conditions)))
})

test_that("the counts table is random and reproducible", {
  expect_identical(simulateRnaSeqData(), simulateRnaSeqData())
  expect_different(simulateRnaSeqData(seed = 0L), simulateRnaSeqData(seed = 1L))
})

test_that("calling the function doesn't change the global state of the random number generator", {
  expect_rng_unchanged(simulateRnaSeqData)
})
