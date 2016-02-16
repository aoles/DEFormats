library("DESeq2")

context("Coerce objects to DGEList")

se = mockRnaSeqData(output = "RangedSummarizedExperiment")
dds = DESeqDataSet(se, design = ~ condition)

test_that("result is of type DGEList", {
  expect_is(as.DGEList(dds), "DGEList")
  expect_is(as(dds, "DGEList"), "DGEList")
})
