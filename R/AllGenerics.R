#' Convert to DGEList
#' 
#' Coerces an object to \linkS4class{DGEList}.
#' 
#' @param x an \R object
#' @param ... additional arguments to be passed to methods
#' @return A \code{\linkS4class{DGEList}} object.
#' @template author
#' @export
as.DGEList = function (x, ...) UseMethod("as.DGEList")
#setGeneric ("as.DGEList", function (x, ...) standardGeneric("as.DGEList") )

#' Convert to DESeqDataSet
#' 
#' Coerces an object to \linkS4class{DESeqDataSet}.
#' 
#' @inheritParams as.DGEList
#' @return A \code{\linkS4class{DESeqDataSet}} object
#' @template author
#' @export
as.DESeqDataSet = function (x, ...) UseMethod("as.DESeqDataSet")
#setGeneric ("as.DESeqDataSet", function (x, ...) standardGeneric("as.DESeqDataSet")

#' DGEList Constructor Generic
#'
#' Creates a \linkS4class{DGEList} object.
#' 
#' @template author
#' @export
setGeneric("DGEList", valueClass = "DGEList")
