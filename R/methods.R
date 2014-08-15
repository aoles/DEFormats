## DESeqDataSet -> DGEList

setAs(from = "DESeqDataSet",
      to = "DGEList",
      def = function(from) {
        ## get annotation information from GRangesList
        genes = as.data.frame(rowData(from))
        if ( nrow(genes) == 0 ) genes = NULL
        
        to = DGEList(
          counts = counts(from), 
          group = from$group,
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
#' data(mockRnaSeqData)
#' group = rep(c("case", "control"), each = 3)
#' 
#' dds = DESeqDataSetFromMatrix(
#'   countData = mockRnaSeqData,
#'   colData = as.data.frame(group, row.names = colnames(mockRnaSeqData)),
#'   design = ~ group
#'   )
#'   
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
#' data(mockRnaSeqData)
#' group = rep(c("case", "control"), each = 3)
#' 
#' dgelist <- DGEList(counts = mockRnaSeqData, group = group)
#' dgelist
#' 
#' as.DESeqDataSet(dgelist)
setMethod ("as.DESeqDataSet", signature(x = "DGEList"), function(x) as(x, "DESeqDataSet") )
