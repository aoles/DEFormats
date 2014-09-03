#' Convert between differential gene expression data formats
#' 
#' \pkg{DEFormats} facilitates data conversion between formats used by gene 
#' expression analysis packages.
#' 
#' Currently the package supports data conversion between \pkg{DESeq2} and 
#' \pkg{edgeR}, i.e. between the \code{\link[DESeq2]{DESeqDataSet-class}} and
#' \code{\link[edgeR]{DGEList-class}}, objects, respectively.
#' 
#' Objects can be coerced using the following methods 
#' \itemize{
#'   \item \code{\link{as.DESeqDataSet}}
#'   \item \code{\link{as.DGEList}}
#' }
#' @author Andrzej Oles \email{andrzej.oles@@embl.de}
#' @docType package
#' @name DEFormats
NULL
