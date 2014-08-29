#' Example counts table of RNA-seq data
#' 
#' Simulated expression data of an RNA-seq experiment.
#' 
#' The count table is generated using the \code{\link[DESeq2]{makeExampleDESeqDataSet}} method from the \pkg{DESeq2} package.
#' 
#' @param output Output format
#' @param rep Number of replicates for each condition
#' @param conditions Condition names
#' @param ... Arguments passed to \code{\link[DESeq2]{makeExampleDESeqDataSet}}
#' @return Depending on the \code{output} setting a matrix or an \code{\link[GenomicRanges]{SummarizedExperiment-class}} object.
#' @author Andrzej Oles <\email{andrzej.oles@@embl.de}>
#' @examples
#' ## count data
#' mockRnaSeqData = mockRnaSeqData()
#' 
#' ## return an SummarizedExperiment object
#' mockRnaSeqDataSE = mockRnaSeqData(output = "SummarizedExperiment")
#' 
#' require("GenomicRanges")
#' identical(mockRnaSeqData, assay(mockRnaSeqDataSE))
mockRnaSeqData = function(output = c("matrix", "SummarizedExperiment"), rep = 3, conditions = c("Case", "Control"), ...) {
  output = match.arg(output)
  
  set.seed(1)
  dds = makeExampleDESeqDataSet(m = 2 * rep, ...)
  colnames(dds) = paste(rep(conditions, each = rep), rep(seq_len(rep), 2))
  
  switch (output, 
          matrix = counts(dds), 
          SummarizedExperiment = as(dds, "SummarizedExperiment")
  )  
}

# save(.mockRnaSeqData(), file = "./data/mockRnaSeqData.RData", compress = "xz")
