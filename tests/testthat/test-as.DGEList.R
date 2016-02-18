library("DESeq2")

context("Coerce objects to DGEList")

se = mockRnaSeqData(output = "RangedSummarizedExperiment")
dds = DESeqDataSet(se, design = ~ condition)
dge = as.DGEList(dds)

test_that("result is of type DGEList", {
  expect_is(dge, "DGEList")
  expect_is(as(dds, "DGEList"), "DGEList")
})

test_that("the counts table is preserved", {
  expect_identical(counts(dds), dge$counts)
})
