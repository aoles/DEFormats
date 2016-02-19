library("DESeq2")

context("Coerce objects to DESeqDataSet")

## construct DGEList object
se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
m = ncol(se)
dge = DGEList(assay(se),
        lib.size = rep(40000L, length = ncol(se)),
        norm.factors = rep(1:2, length = ncol(se)),
        group = factor(rep(c("a", "b"), length = ncol(se))),
        genes = as.data.frame(rowRanges(se)))
        
## add offsets
dge$offset = log(simulateNormFactors(dimnames = dimnames(dge$counts)))

## coerce to DGEList
dds = as.DESeqDataSet(dge)

test_that("result is of type DESeqDataSet", {
  expect_is(dds, "DESeqDataSet")
})

test_that("elements from DGEList are carried over to DESeqDataSet", {
  expect_identical(dge$counts, counts(dds))
  expect_identical(dge$samples$group, colData(dds)$group)
  expect_identical(dge$samples$condition, colData(dds)$condition)
  expect_identical(dge$samples$norm.factors, colData(dds)$norm.factors)
  expect_identical(dge$samples$lib.size, colData(dds)$lib.size)
  expect_identical(dge$genes, as.data.frame(rowRanges(dds)))
})

test_that("converting to DESeqDataSet and back to DGEList restores the original object", {
  dgedds = as.DGEList(dds)
  expect_identical(dge$counts, dgedds$counts)
  expect_identical(dge$samples, dgedds$samples)
  expect_identical(dge$genes, dgedds$genes)
  expect_identical(dge$offset, dgedds$offset)
})
