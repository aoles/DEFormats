library("DESeq2")

context("DGEList Constructor Generic")

se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
mx = assay(se)

test_that("the returned value is of class DGEList", {
  expect_is(DGEList(), "DGEList")
  expect_is(DGEList(se), "DGEList")
  expect_is(DGEList(mx), "DGEList")
})

test_that("the resulting objects contain the same information", {
  condition = colData(se)$condition
  genes = as.data.frame(rowRanges(se))
  expect_identical(DGEList(se, group = condition),
                   DGEList(mx, group = condition, genes = genes))
})
