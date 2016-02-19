#' Simulate Normalization Factors
#' 
#' Simulate gene-specific normalization factors for each sample of an RNA-seq experiment.
#' 
#' @inheritParams simulateRnaSeqData
#' @param ... arguments passed to \code{\link[base]{matrix}}
#' @return A matrix with \code{n} rows and \code{m} columns containing the normalization factors.
#' @seealso simulateRnaSeqData
#' @template author
#' @examples
#' require("DESeq2")
#' 
#' ## normalization factors
#' se = simulateRnaSeqData(output = "RangedSummarizedExperiment")
#' 
#' dds = DESeqDataSet(se, design = ~ condition)
#' 
#' normalizationFactors(dds) = simulateNormFactors()
#' @export
simulateNormFactors = function(n = 1000L, m = 6L, seed = 0L, ...) {
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
  
  normFactors = matrix(runif(n*m, 0.5, 1.5), n, m, ...)
  
  # the normalization factors matrix should not have 0's in it
  # it should have geometric mean near 1 for each row
  normFactors = normFactors / exp(rowMeans(log(normFactors)))
  
  normFactors
}
