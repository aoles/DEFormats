## DESeqDataSet -> DGEList

setAs(from = "DESeqDataSet",
      to = "DGEList",
      def = function(from) {
        ## get annotation information from GRangesList
        genes = as.data.frame(rowData(from))
        if ( nrow(genes) == 0 ) genes = NULL
        
        ## group by last variable
        v = as.character(attr(terms(design(from)), "variables"))
        group = v[length(v)]
        
        to = DGEList(
          counts = counts(from),
          group = colData(from)[[group]],   
          genes = genes
          )
                
        ## copy normalization factors if present
        if ( !is.null( (nf = normalizationFactors(from)) ) )
          to$offset = log(nf)
        
        return(to)
      })

#' @describeIn as.DGEList Coerce \code{\link[DESeq2]{DESeqDataSet-class}} objects to \code{\link[edgeR]{DGEList-class}}.
#' @examples
#' require("DESeq2")
#' 
#' se = mockRnaSeqData(output = "SummarizedExperiment")
#' se
#' 
#' dds = DESeqDataSet(se, design = ~ condition)   
#' dds
#' 
#' as.DGEList(dds)
setMethod ("as.DGEList", signature(x = "DESeqDataSet"), function(x) as(x, "DGEList") )


## DGEList -> DESeqDataSet

setAs(from = "DGEList",
      to = "DESeqDataSet",      
      def = function(from) {
        to = DESeqDataSetFromMatrix(
          countData = from$counts,
          colData = from$samples["group"],
          design = formula("~ group", env = globalenv())
          )
        
        ## copy normalization factors if present
        if ( !is.null( (nf = from$offset) ) )
          normalizationFactors(to) = exp(nf)
        
        return(to)
      })

#' @describeIn as.DESeqDataSet Coerce \code{\link[edgeR]{DGEList-class}} objects to \code{\link[DESeq2]{DESeqDataSet-class}}.
#' @examples
#' require("edgeR")
#' 
#' counts = mockRnaSeqData()
#' group = rep(c("case", "control"), each = 3)
#' 
#' dge = DGEList(counts = counts, group = group)
#' dge
#' 
#' as.DESeqDataSet(dge)
setMethod ("as.DESeqDataSet", signature(x = "DGEList"), function(x) as(x, "DESeqDataSet") )


#' DGEList Constructor Generic
#' 
#' Create a \code{\link[edgeR]{DGEList-class}} object 
#' @rdname DGEList
setMethod ("DGEList", signature(counts = "SummarizedExperiment"), function(
  counts = SummarizedExperiment(),
  lib.size = colSums(assay(counts)),
  norm.factors = rep(1, ncol(counts)),
  group = rep(1, ncol(counts)),
  genes = as.data.frame(rowData(counts)),
  remove.zeros = FALSE
  ) DGEList(
    assay(counts),
    lib.size = lib.size,
    norm.factors = norm.factors,
    group = group,
    genes = genes,
    remove.zeros = remove.zeros
    ) 
)
