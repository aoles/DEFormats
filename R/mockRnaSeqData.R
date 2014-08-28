#' Example counts table of RNA-seq data
#' 
#' Simulated expression data of an RNA-seq experiment.
#' 
#' @usage data("mockRnaSeqData")
#' @examples
#' data("mockRnaSeqData")
#' str(mockRnaSeqData)
#' @format A 1000 x 6 integer matrix with rows representing genes and columns different samples.
#' @source The count table was generated using the \code{\link[DESeq2]{makeExampleDESeqDataSet}} method from the \pkg{DESeq2} package.
#' @name mockRnaSeqData
NULL

.mockRnaSeqData = function(rep = 3) {
  set.seed(1)
  dds = makeExampleDESeqDataSet(m = 2 * rep)
  colnames(dds) = paste(rep(c("Case", "Control"), each = rep), rep(seq_len(rep), 2))
  counts(dds)
}

# save(.mockRnaSeqData(), file = "./data/mockRnaSeqData.RData", compress = "xz")
