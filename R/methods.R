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

## set method
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

## set method
setMethod ("as.DESeqDataSet", signature(x = "DGEList"), function(x) as(x, "DESeqDataSet") )
