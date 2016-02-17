library("DESeq2")

context("DGEList Constructor Generic")

se = mockRnaSeqData(output = "RangedSummarizedExperiment")
mx = assay(se)

test_that("the returned value is of class DGEList", {
  expect_is(DGEList(), "DGEList")
  expect_is(DGEList(se), "DGEList")
  expect_is(DGEList(mx), "DGEList")
})

test_that("the resulting objects contain the same information", {
  colData(se)$condition = NULL
  expect_identical(DGEList(se, genes = NULL), DGEList(mx))
})
