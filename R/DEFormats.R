#' Convert between differential gene expression data formats
#' 
#' \pkg{DEFormats} facilitates data conversion between various formats used by
#' different gene expression analysis packages.
#' 
#' Currently the package supports data conversion between \pkg{DESeq2} and 
#' \pkg{edgeR}, i.e. between the \code{\link[DESeq2]{DESeqDataSet-class}} and 
#' \code{\link[edgeR]{DGEList-class}}, objects, respectively.
#' 
#' Objects can be coerced using the following methods \itemize{ \item
#' \code{\link{as.DESeqDataSet}} \item \code{\link{as.DGEList}} }
#' @author Andrzej Oleś \email{andrzej.oles@@embl.de}
#' @docType package
#' @name DEFormats
#' @importFrom edgeR DGEList
#' @importClassesFrom edgeR DGEList
#' @importFrom DESeq2 DESeqDataSetFromMatrix normalizationFactors counts design normalizationFactors<- makeExampleDESeqDataSet
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importFrom SummarizedExperiment rowRanges colData assay as.data.frame SummarizedExperiment
#' @importFrom methods setAs setGeneric setMethod
NULL
