## does not need to be exported
setAs(from = "DESeqDataSet",
      to = "DGEList",
      def = function(from) {
        
        ## get annotation information from GRangesList
        genes = as.data.frame(rowRanges(from))
        if ( nrow(genes) == 0 ) genes = NULL
        
        ## get sample description
        samples = as.data.frame(colData(from))
        
        ## group by last variable
        v = as.character(attr(terms(design(from)), "variables"))
        group = v[length(v)]
        
        ## count matrix
        counts = counts(from)
        
        to = DGEList(
          counts = counts,
          lib.size = if ( is.null( (ls = samples$lib.size) ) ) colSums(counts) else ls,
          norm.factors = if ( is.null( (nf = samples$norm.factors) ) ) rep(1, ncol(counts)) else nf,
          group = samples[[group]],
          genes = genes
          )
        
        ## copy remaining sample metadata removing any duplicates
        df = data.frame(to$samples, samples)    
        to$samples = df[!duplicated(as.list(df))]
        
        ## copy normalization factors if present
        if ( !is.null( (nf = normalizationFactors(from)) ) )
          to$offset = log(nf)
        
        return(to)
      })

#' @describeIn as.DGEList Coerce \code{\linkS4class{DESeqDataSet}} objects to \code{\link[edgeR]{DGEList-class}}.
#' @examples
#' require("DESeq2")
#' 
#' se = mockRnaSeqData(output = "RangedSummarizedExperiment")
#' se
#' 
#' dds = DESeqDataSet(se, design = ~ condition)   
#' dds
#' 
#' as.DGEList(dds)
#' @export
as.DGEList.DESeqDataSet = function (x, ...) as(x, "DGEList")

## does not need to be exported
setAs(from = "DGEList",
      to = "DESeqDataSet",      
      def = function(from) {
        to = DESeqDataSetFromMatrix(
          countData = from$counts,
          colData = from$samples,
          design = formula("~ group", env = globalenv())
          )
        
        ## copy normalization factors if present
        if ( !is.null( (nf = from$offset) ) )
          normalizationFactors(to) = exp(nf)
        
        return(to)
      })

#' @describeIn as.DESeqDataSet Coerce \code{\link[edgeR]{DGEList-class}} objects to \code{\linkS4class{DESeqDataSet}}.
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
#' @export
as.DESeqDataSet.DGEList = function (x, ...) as(x, "DESeqDataSet")

#' @param counts read counts, either a numeric matrix or a \linkS4class{RangedSummarizedExperiment} object
#' @param lib.size numeric vector giving the total count (sequence depth) for each library
#' @param norm.factors numeric vector of normalization factors that modify the library sizes
#' @param group vector or factor giving the experimental group/condition for each sample/library
#' @param genes data frame containing annotation information for each gene
#' @param remove.zeros logical, whether to remove rows that have 0 total count 
#' 
#' @rdname DGEList
#' @export
setMethod ("DGEList", signature(counts = "RangedSummarizedExperiment"),
  function(counts,
           lib.size = colSums(assay(counts)),
           norm.factors = rep(1, ncol(counts)),
           group = rep(1, ncol(counts)),
           genes = as.data.frame(rowRanges(counts)),
           remove.zeros = FALSE) {
    dge = DGEList(assay(counts),
                  lib.size = lib.size,
                  norm.factors = norm.factors,
                  group = group,
                  genes = genes,
                  remove.zeros = remove.zeros)
    
    ## copy remaining sample metadata removing any duplicates
    df = data.frame(dge$samples, colData(counts))
    dge$samples = df[!duplicated(as.list(df))]
    
    dge
  }
)
