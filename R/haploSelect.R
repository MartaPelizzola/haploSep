#' Model selector for haplotype separation
#'
#' More details
#'
#' @param data Numeric observation matrix with nrow(data) being the number of
#' SNP locations and ncol(data) the number of samples (e.g. time points).
#'
#' @return integer. The estimated number of essential haplotypes.
#'
#' @examples
#'  data(ExampleDataset)
#'  haploSelect(data = Y)
#' @export

haploSelect = function(data) {
  M = ncol(data)
  n = nrow(data)
  
  # default parameters
  m.min = 2          # Integer value which gives minimal number of haplotypes
  m.max = min(M, 10) # Integer value which gives maximal number of haplotypes

  # SVD
  sv = svd(data)$d
  # Gavish and Donoho criterion
  beta = M/n
  w    = 0.56 * beta^3 - 0.95 * beta^2 + 1.82 * beta + 1.43
  tau  = w * median(sv)
  ans  = max(sum(sv >= tau),2)
  if (ans > m.max) {
    return(m.max)
  } else if (ans < m.min) {
    return(m.min)
  } else {
    return(ans)
  }
}
