library("DESeq2")

context("Coerce objects to DGEList")

rep = 3L
se = mockRnaSeqData(output = "RangedSummarizedExperiment", rep = rep)
dds = DESeqDataSet(se, design = ~ condition)
DF = colData(dds)
DF$norm.factors = 1:2*rep
DF$lib.size = rep(40000L, 2*rep)
DF$group = rep(c("a", "b"), each = rep)
colData(dds) = DF

dge = as.DGEList(dds)

test_that("result is of type DGEList", {
  expect_is(dge, "DGEList")
  expect_is(as(dds, "DGEList"), "DGEList")
})

test_that("elements from DESeqDataSet are carried over to DGEList", {
  expect_identical(dge$counts, counts(dds))
  expect_identical(dge$samples$group, colData(dds)$condition)
  expect_identical(dge$samples$group.1, colData(dds)$group)
  expect_identical(dge$samples$norm.factors, colData(dds)$norm.factors)
  expect_identical(dge$samples$lib.size, colData(dds)$lib.size)
  expect_identical(dge$genes, as.data.frame(rowRanges(dds)))
})
