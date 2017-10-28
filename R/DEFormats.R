#' Convert Between Differential Gene Expression Data Formats
#' 
#' \pkg{DEFormats} provides data converters between various formats used by 
#' different gene expression analysis packages.
#' 
#' Currently the package supports data conversion between \pkg{DESeq2} and 
#' \pkg{edgeR}, i.e., between \code{\linkS4class{DESeqDataSet}} and 
#' \code{\linkS4class{DGEList}} objects, respectively.
#' 
#' Objects can be coerced using the following methods \itemize{ \item 
#' \code{\link{as.DESeqDataSet}} \item \code{\link{as.DGEList}} }
#' @template author
#' @docType package
#' @name DEFormats
#' @import data.table
#' @importFrom checkmate qassert
#' @importFrom edgeR DGEList
#' @importClassesFrom edgeR DGEList
#' @importFrom DESeq2 DESeqDataSetFromMatrix normalizationFactors counts design
#'   normalizationFactors<- makeExampleDESeqDataSet
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importFrom SummarizedExperiment rowRanges colData assay as.data.frame
#'   SummarizedExperiment
#' @importFrom GenomicRanges GRanges makeGRangesListFromDataFrame
#' @importFrom methods as is new setAs setGeneric setMethod
#' @importFrom stats formula runif terms
NULL
