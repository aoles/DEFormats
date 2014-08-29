## set S4 generics

#' Coerce objects to DGEList
#' 
#' Coerce objects to \code{\link[edgeR]{DGEList-class}}.
#' 
#' @param x An object
#' @param ... Additional arguments to be passed to methods
#' @return A \code{\link[edgeR]{DGEList-class}} object
#' @author Andrzej Oles \email{andrzej.oles@@embl.de}
setGeneric ("as.DGEList", function (x, ...) standardGeneric("as.DGEList") )

#' Coerce objects to DESeqDataSet
#' 
#' Coerce objects to \code{\link[DESeq2]{DESeqDataSet-class}}.
#' 
#' @inheritParams as.DGEList
#' @return A \code{\link[DESeq2]{DESeqDataSet-class}} object
#' @author Andrzej Oles \email{andrzej.oles@@embl.de}
setGeneric ("as.DESeqDataSet", function (x, ...) standardGeneric("as.DESeqDataSet") )

# #' DGEList Constructor Generic
# #' 
# #' Create a \code{\link[edgeR]{DGEList-class}} object
#setGeneric("DGEList")
