library("DESeq2")

context("Coerce objects to DGEList")

rep = 3L
se = mockRnaSeqData(output = "RangedSummarizedExperiment", rep = rep)
DF = colData(se)
names(DF) = "group"
DF$lib.size = 40000L
DF$norm.factors = 1:2
DF$condition = rep(factor(c("a", "b")), rep)
colData(se) = DF

dds = DESeqDataSet(se, design = ~ condition + group)
dge = as.DGEList(dds)
dgedds = as.DESeqDataSet(dge)

test_that("result is of type DGEList", {
  expect_is(dge, "DGEList")
  expect_is(as(dds, "DGEList"), "DGEList")
})

test_that("elements from DESeqDataSet are carried over to DGEList", {
  expect_identical(dge$counts, counts(dds))
  expect_identical(dge$samples$group, colData(dds)$group)
  expect_identical(dge$samples$condition, colData(dds)$condition)
  expect_identical(dge$samples$norm.factors, colData(dds)$norm.factors)
  expect_identical(dge$samples$lib.size, colData(dds)$lib.size)
  expect_identical(dge$genes, as.data.frame(rowRanges(dds)))
})

test_that("converting to DGEList and back to DESeqDataSet restores the original object", {
  expect_identical(counts(dds), counts(dgedds))
  expect_identical(colData(dds), colData(dgedds))
  expect_identical(rowRanges(dds), rowRanges(dgedds))
  #expect_identical(normalizationFactors(dge), normalizationFactors(dgedds))
})
