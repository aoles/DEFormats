## does not need to be exported
setAs(from = "DESeqDataSet",
      to = "DGEList",
      def = function(from) {
        
        ## argument list to DGEList constructor
        args = list()
        
        ## count matrix
        args$counts = counts(from)
        
        ## get annotation information from GRangesList
        genes = as.data.frame(rowRanges(from))
        if ( nrow(genes) == nrow(args$counts) ) 
          args$genes = genes
        
        ## get sample description
        samples = as.data.frame(colData(from))
        colnames = names(samples)
        
        add_arg = function(col, arg=col) {
          if (col %in% colnames) {
            args[[arg]] <<- samples[[col]]
            colnames <<- colnames[colnames!=col]
          }
        }
        
        add_arg("lib.size")
        
        add_arg("norm.factors")
        
        ## group by last variable
        v = as.character(attr(terms(design(from)), "variables"))
        group = v[length(v)]
        
        add_arg(group, "group")
        
        ## any remaining sample metadata
        args$samples = samples[colnames]
        
        to = do.call("DGEList", args)
        
        ## copy normalization factors if present
        if ( !is.null( (nf = normalizationFactors(from)) ) )
          to$offset = log(nf)
        
        return(to)
      })

#' @describeIn as.DGEList Coerce \code{\linkS4class{DESeqDataSet}} objects to
#'   \code{\link[edgeR]{DGEList-class}}.
#' @examples
#' require("DESeq2")
#' 
#' se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
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
        args = list()
        
        args$countData = from$counts
        args$colData = from$samples
        args$design = design = formula("~ group", env = .GlobalEnv)
        
        if (!is.null((genes = from$genes))) {
          ## attempt first constructing a GenomicRanges object
          genes = tryCatch(
            GRanges(genes),
            ## if this fails try creating a GenomicRangesList from feature fragments
            error = function(e)
              tryCatch({
                ## use data.table to unpack feature fragments
                ## https://stackoverflow.com/a/43431847/2792099
                gdt = as.data.table(genes, keep.rownames = "ID")
                sep = getOption("DEFormats.featureFragments.sep", ",")
                gdt = gdt[, lapply(.SD, function(x)
                  unlist(strsplit(as.character(x), sep, fixed = TRUE))), by = "ID"]
                makeGRangesListFromDataFrame(gdt, split.field = "ID")
              },
              ## coercion to GRanges or GRangesList failed, attach as is
              error = function(e) genes ) )
          if (is(genes, "GenomicRanges_OR_GRangesList"))
            args$rowRanges = genes
          else
            args$rowData = genes
        }
        
        to = do.call("DESeqDataSetFromMatrix", args)
        
        ## copy normalization factors if present
        if (!is.null( (nf = from$offset) ))
          normalizationFactors(to) = exp(nf)
        
        return(to)
      })

#' @describeIn as.DESeqDataSet Coerce \code{\link[edgeR]{DGEList-class}} objects
#'   to \code{\linkS4class{DESeqDataSet}}.
#' @examples
#' require("edgeR")
#' 
#' counts = simulateRnaSeqData()
#' group = rep(c("case", "control"), each = 3)
#' 
#' dge = DGEList(counts = counts, group = group)
#' dge
#' 
#' as.DESeqDataSet(dge)
#' @export
as.DESeqDataSet.DGEList = function (x, ...) as(x, "DESeqDataSet")

#' @param counts read counts, either a numeric matrix or a
#'   \linkS4class{RangedSummarizedExperiment} object
#' @param lib.size numeric vector giving the total count (sequence depth) for
#'   each library
#' @param norm.factors numeric vector of normalization factors that modify the
#'   library sizes
#' @param samples data frame containing information for each sample
#' @param group vector or factor giving the experimental group/condition for
#'   each sample/library
#' @param genes data frame containing annotation information for each gene
#' @param remove.zeros logical, whether to remove rows that have 0 total count
#' @examples 
#' se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
#' 
#' ## Initialize a DGEList from a RangedSummarizedExperiment object
#' DGEList(se)
#' @rdname DGEList
#' @export
setMethod ("DGEList", signature(counts = "RangedSummarizedExperiment"),
  function(counts = new("RangedSummarizedExperiment"),
           lib.size = colData(counts)$lib.size,
           norm.factors = colData(counts)$norm.factors,
           samples = colData(counts),
           group = colData(counts)$group,
           genes = as.data.frame(rowRanges(counts)),
           remove.zeros = FALSE) {
      # remove duplicated columns
      if (!is.null(samples)) {
        dup = names(samples) %in% c("lib.size", "norm.factors", "group")
        samples = if (all(dup)) {
          NULL
        } else {
          samples = samples[!dup]
        }
      }
      DGEList(assay(counts),
              lib.size = lib.size,
              norm.factors = norm.factors,
              samples = samples,
              group = group,
              genes = genes,
              remove.zeros = remove.zeros)
  }
)
