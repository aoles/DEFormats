#' Example counts table of RNA-seq data
#' 
#' Simulated expression data of an RNA-seq experiment.
#' 
#' The count table is generated using the
#' \code{\link[DESeq2]{makeExampleDESeqDataSet}} method from the \pkg{DESeq2}
#' package.
#' 
#' @param output output type
#' @param n number of genes
#' @param m number of samples
#' @param seed a single integer value specifying the random number generator
#'   seed
#' @param ... arguments passed to \code{\link[DESeq2]{makeExampleDESeqDataSet}}
#' @return Depending on the \code{output} setting a matrix or an
#'   \code{\linkS4class{RangedSummarizedExperiment}} object.
#' @seealso simulateNormFactors
#' @template author
#' @examples
#' ## count data matrix
#' mx = simulateRnaSeqData()
#' head(mx)
#' 
#' ## return an RangedSummarizedExperiment object
#' se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
#' se
#' @export
simulateRnaSeqData = function(output = c("matrix","RangedSummarizedExperiment"),
                              n = 1000, m = 6, seed = 0L, ...) {
  output = match.arg(output)
  qassert(n, "X1[1,)")
  qassert(m, "X1[1,)")
  qassert(seed, "X1")
  
  oldseed <- .GlobalEnv$.Random.seed
  on.exit({
    if (!is.null(oldseed)) 
      .GlobalEnv$.Random.seed <- oldseed
    else
      rm(".Random.seed", envir = .GlobalEnv)
    })
  set.seed(seed)
  
  dds = makeExampleDESeqDataSet(n = n, m = m, ...)
  
  switch (output, 
          matrix = counts(dds), 
          RangedSummarizedExperiment = as(dds, "RangedSummarizedExperiment")
  )  
}
