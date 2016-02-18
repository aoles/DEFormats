#' Simulate Normalization Factors
#' 
#' Simulate gene-specific normalization factors for each sample of an RNA-seq experiment.
#' 
#' @param n number of genes
#' @param m number of samples
#' @param seed random number generator seed
#' @return A \code{n} times \code{m} matrix containing the normalization factors.
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
simulateNormFactors = function(n = 1000L, m = 6L, seed = 0L) {
  oldseed <- .GlobalEnv$.Random.seed
  on.exit({
    if (!is.null(oldseed)) 
      .GlobalEnv$.Random.seed <- oldseed
    else
      rm(".Random.seed", envir = .GlobalEnv)
    })
  set.seed(seed)
  
  normFactors = matrix(runif(n*m, 0.5, 1.5), n, m)
  
  # the normalization factors matrix should not have 0's in it
  # it should have geometric mean near 1 for each row
  normFactors = normFactors / exp(rowMeans(log(normFactors)))
  
  normFactors
}
