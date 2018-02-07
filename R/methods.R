.featureFragments.sep = ","

## does not need to be exported
setAs(from = "DESeqDataSet",
      to = "DGEList",
      def = function(from) {
        
        ## argument list to DGEList constructor
        args = list()
        
        ## count matrix
        args$counts = counts(from)
        
        ## get row metadata
        genes = as.data.frame(rowRanges(from))
        rownames = NULL
        
        if (nrow(genes) > 1L) {
          if (nrow(genes) > nrow(args$counts)) {
            ## pack GRangesList into csv-style data.frame
            sep = getOption("DEFormats.featureFragments.sep", .featureFragments.sep)
            setDT(genes)
            genes = genes[, lapply(.SD, function(x) paste(x, collapse = sep)),
                          .SDcols = -c("group", "group_name"),
                          by = "group_name"]
            rownames = "group_name"
          }
        }
        else {
          genes = rowData(from)
          if (!is.null(genes$rownames))
            rownames = "rownames"
        }
        
        ## attach non-empty features metadata
        if (length(genes) > 0L)
          args$genes = data.frame(genes, row.names = rownames)
        
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
        
        has_mcols = FALSE
        
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
                sep = getOption("DEFormats.featureFragments.sep", .featureFragments.sep)
                gdt = gdt[, lapply(.SD, function(x)
                  unlist(strsplit(as.character(x), sep, fixed = TRUE))), by = "ID"]
                makeGRangesListFromDataFrame(gdt, split.field = "ID")
              },
              ## coercion to GRanges or GRangesList failed, attach as DataFrame
              error = function(e) {
                ## add row names column as they are not supported in rowData
                rownames = attr(genes, "row.names")
                if (identical(rownames, (seq_len(nrow(genes)))))
                  DataFrame(genes)
                else
                  DataFrame(rownames, genes)
              } ) )
          if (is(genes, "GenomicRanges_OR_GRangesList"))
            args$rowRanges = genes
          else
            has_mcols = TRUE
        }
        
        to = do.call("DESeqDataSetFromMatrix", args)
        
        ## attach rows metadata which could not be coerced to GenomicRanges or
        ## GRangesList as adviced in ?DESeqDataSet
        if (has_mcols)
          mcols(to) = genes
        
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

#' @inheritParams edgeR::DGEList
#' @param counts read counts, either a numeric matrix or a
#'   \linkS4class{RangedSummarizedExperiment} object.
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
           group = NULL,
           genes = as.data.frame(rowRanges(counts)),
           remove.zeros = FALSE) {
      # remove duplicated columns
      if (!is.null(samples)) {
        dup = names(samples) %in% c("lib.size", "norm.factors")
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
