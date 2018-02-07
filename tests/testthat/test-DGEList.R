library("DESeq2")

context("DGEList Constructor Generic")

se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
mx = assay(se)
names(colData(se)) = "group"

test_that("the returned value is of class DGEList", {
  expect_is(DGEList(), "DGEList")
  expect_is(DGEList(se), "DGEList")
  expect_is(DGEList(mx), "DGEList")
})

test_that("the resulting objects contain the same information", {
  samples = colData(se)
  genes = as.data.frame(rowRanges(se))
  dge = DGEList(se)
  expect_identical(dge, DGEList(mx, samples = samples, genes = genes))
  rse = as(as.DESeqDataSet(dge), "RangedSummarizedExperiment")
  expect_identical(dge, DGEList(rse))
})
