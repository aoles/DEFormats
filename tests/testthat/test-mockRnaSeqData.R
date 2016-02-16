library("DESeq2")

context("Simulate RNA-seq data")

se = mockRnaSeqData(output = "RangedSummarizedExperiment")
dds = DESeqDataSet(se, design = ~ condition)

test_that("output is of desired type", {
  expect_is(mockRnaSeqData(output = "matrix"), "matrix")
  expect_is(mockRnaSeqData(output = "RangedSummarizedExperiment"), "RangedSummarizedExperiment")
})

test_that("the number of replicates is correct", {
  rep = 5L
  expect_equal(ncol(mockRnaSeqData(rep = rep)), 2*rep)
})

test_that("columns are named after conditions", {
  conditions = c("a", "b")
  n = colnames(mockRnaSeqData(conditions = conditions))
  n = sub(" [0-9]+$", "",  n)
  expect_identical(n , rep(conditions, each = length(n)/length(conditions)))
})

test_that("the counts table is random and reproducible", {
  expect_identical(mockRnaSeqData(), mockRnaSeqData())
  expect_false(identical(mockRnaSeqData(seed = 0L), mockRnaSeqData(seed = 1L)))
})

test_that("calling the function doesn't change the state of random number generator", {
  oldseed <- .GlobalEnv$.Random.seed
  mockRnaSeqData()
  expect_identical(oldseed, .GlobalEnv$.Random.seed)
})
