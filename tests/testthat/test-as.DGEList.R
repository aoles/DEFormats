library("DESeq2")

context("Coerce objects to DGEList")

## construct DESeqDataSet object
se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
DF = colData(se)
DF$group = factor(c("a", "b"))
DF$lib.size = 40000L
DF$norm.factors = 1:2
colData(se) = DF[c(2L:ncol(DF), 1L)] # swap columns to ensure reproducibility with DGEList 
dds = DESeqDataSet(se, design = ~ condition + group)

## add normalization Factors
normalizationFactors(dds) = simulateNormFactors()

## coerce to DGEList
dge = as.DGEList(dds)

test_that("result is of type DGEList", {
  expect_is(dge, "DGEList")
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
  dgedds = as.DESeqDataSet(dge)
  expect_identical(counts(dds), counts(dgedds))
  expect_identical(colData(dds), colData(dgedds))
  expect_identical(rowRanges(dds), rowRanges(dgedds))
  expect_identical(normalizationFactors(dds), normalizationFactors(dgedds))
})
