#' Example counts table of RNA-seq data
#' 
#' Simulated expression data of an RNA-seq experiment.
#' 
#' The count table is generated using the \code{\link[DESeq2]{makeExampleDESeqDataSet}} method from the \pkg{DESeq2} package.
#' 
#' @param output Output format
#' @param rep Number of replicates for each condition
#' @param conditions Condition names
#' @param seed random number generator seed
#' @param ... Arguments passed to \code{\link[DESeq2]{makeExampleDESeqDataSet}}
#' @return Depending on the \code{output} setting a matrix or an \code{\link[SummarizedExperiment]{RangedSummarizedExperiment-class}} object.
#' @template author
#' @examples
#' ## count data
#' mockRnaSeqData = mockRnaSeqData()
#' 
#' ## return an RangedSummarizedExperiment object
#' mockRnaSeqDataSE = mockRnaSeqData(output = "RangedSummarizedExperiment")
#' 
#' require("SummarizedExperiment")
#' identical(mockRnaSeqData, assay(mockRnaSeqDataSE))
#' @export
mockRnaSeqData = function(output = c("matrix", "RangedSummarizedExperiment"), rep = 3, conditions = c("Case", "Control"), seed = 0L, ...) {
  output = match.arg(output)
  
  set.seed(seed)
  
  dds = makeExampleDESeqDataSet(m = 2 * rep, ...)
  colnames(dds) = paste(rep(conditions, each = rep), rep(seq_len(rep), 2))
  
  switch (output, 
          matrix = counts(dds), 
          RangedSummarizedExperiment = as(dds, "RangedSummarizedExperiment")
  )  
}
