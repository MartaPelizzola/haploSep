#' Separates haplotypes
#'
#' Implementation of interative minimization algorithm for jointly estimating
#' structure of essential haplotypes as well as their relative proportion from
#' allele frequency matrix of mixture.
#'
#' @param data Numeric observation matrix with nrow(data) being the number of
#' SNP locations and ncol(data) the number of samples (e.g. time points).
#' @param nHaplo Integer value which gives the number of essential haplotypes
#' for which haplotype structure and frequency should be estimated from the
#' mixture. When missing, it is estimated via the SVD criterion.
#' @param stabEval logical. If \code{TRUE} (default), the stability analysis
#' of the reconstruction is carried out, otherwise only haplotype structure
#' and frequency are computed.
#' @param bias logical. It indicates whether bias term (constant over SNP
#' locations but different for various samples) should be added. By
#' default, it is \code{TRUE}
#' @param weight matrix of \code{ncol(weight)=ncol(data)} and
#' \code{nrow(weight)=nrow(data)} with weights for the observations least
#' squares optimization step. By default, equal weights are applied.
#' @param nBoot integer. It gives number of bootstrap repetitions for stability
#' score. By default, it takes value \code{15}. It is only needed when
#' \code{stabEval = TRUE}.
#' @param fBoot numeric. Its value lies between 0 and 1 which gives the
#' relative subsample size in the bootstrap runs. By default, it takes
#' value \code{0.95}. It is only needed when \code{stabEval = TRUE}.
#'
#'
#'
#' @return List with entries \code{haploFrq} and \code{haploStr},
#' which is an object from class haplo.
#'
#' \itemize{
#' \item
#' \code{haploFrq} is a matrix with \code{nrow(haploFrq) = nHaplo}
#' and \code{ncol(haploFrq) = ncol(data)} which gives the estimated
#'  frequency of the estimated essential haplotypes.
#' \item
#' \code{haplotStr} is a matrix with \code{nrow(haploStr) = nrow(data)}
#'  and \code{ncol(haploStr) = nHaplo} and entries being either 0 or 1.
#'  It gives the estimated haplotype structure of essential haplotypes.
#' }
#'
#' The returned list has an attribute "\code{nHaplo}", which is the number
#' of essential haplotypes.
#'
#' If \code{stabEval = TRUE} the returned list has three more attributes
#' "\code{R2}", "\code{stabIntFrq}" and "\code{stabScoreStr}".
#'
#' \itemize{
#' \item
#' The attribute "\code{R2}" is an inidicator of how good the model
#' fitting is, in a similar spirit as the \eqn{R^2} for linear models,
#' see \code{\link[stats]{lm}} and \code{\link[stats]{summary.lm}}.
#' \item
#' The attribute "\code{stabIntFrq}" provides a confidence envelope for
#' the estimated haplotype frequency. This is a data frame, containing
#' "\code{lowerBnd}" and "\code{upperBnd}" for each haplotype, which are
#'  0.025 and 0.975 quantiles for the bootstrap samples, respectively.
#' \item
#' The attribute "\code{stabScoreStr}" is the a numeric vector of length
#' \code{nHaplo} with values between 0 and 1.
#' }
#'
#' @examples
#'  1. Reconstruct 5 haplotypes
#'  data(ExampleDataset)
#'  haploSep(data = Y, nHaplo = 5, stabEval = TRUE, bias = TRUE)
#'  
#'  Choose the number of haplotypes to be reconstructed with haploSelect
#'  data(ExampleDataset)
#'  m <- haploSelect(data = Y, bias = TRUE)
#'  haploSep(data = Y, nHaplo = m, stabEval = TRUE, bias = TRUE)
#' @export


haploSep <- function(data, nHaplo, stabEval = TRUE, bias = TRUE, weight = NULL,
                     nBoot = 15, fBoot = 0.95){
  # default parameters & check inputs
  if (!is.matrix(data) || ncol(data) < 2)
    stop("Input 'Data' should be a matrix with at least 2 columns!")
  if (!all(is.finite(data)))
    stop("Input 'data' should contain only finite values!")
  n = nrow(data)
  M = ncol(data)
  if (is.null(weight))
    weight  = matrix(1,n,M)
  al    = c(0,1) # Numeric vector which gives the possible values in the
                 # haplotype structure. By defauls al = c(0,1) corresponding
                 # to allel 0 and 1
  relax = TRUE   # Logical value indicating whether unit sum restriction on
                 # weights should be relaxed.
  tsh   = 1e-3 # Numeric value between 0 and 1 which provides a threshold
               # criterion to the interative Lloyds algorithm. When all
               # weight's increments fall below tsh the algorithm terminates.
               # The default value is tsh=1e-3.
  K     = 3    # Integer value which gives the number of repetitions for the
               # iterative Lloyds algorithm. The final estimate is the one
               # with smallest MSE.
  nBoot = 15   #
  fBoot = 0.95 # Numeric value between 0 and 1 which gives the relative
               # subsample size in the bootstrap runs
  lossFun = "LS" # either "LS" (least squares) or "Likelihood" (binomial)
  z       = NA   # allele counts, needed if lossFun == "Likelihood"

  if (missing(nHaplo))
    nHaplo = haploSelect(data, nHaplo, bias, weight)
  if (nHaplo < 2 || nHaplo > M)
    stop(paste("Input 'nHaplo' should be an integer between 2 and",
               sprintf("# columns (%d) of 'data'!",M)))

  init    = ifelse (length(al)^nHaplo > 250, "random", "combi")
                 # One of "random" or "combi" which determines initialization
                 # of haploFrq. For "combi" the parameter K is silently ignored.

  if (init == "combi") {
    omEst   = estOmega(data, nHaplo, al, bias)
    biasEst = attr(omEst, "bias")
    ans     = lloyd(data, nHaplo, al, tsh, omEst, relax, bias,
                    weight, lossFun, z, biasEst)
  } else {
    ans = vector('list', K)
    for (i in 1:K)
      ans[[i]] = lloyd(data, nHaplo, al, tsh, rOmega(nHaplo,M),
                       relax, bias, weight, lossFun, z)
    se  = sapply(ans, function(x) sum((x$haploStr %*% x$haploFrq + x$bias - data)^2))
    ans = ans[[which.min(se)]]
  }
  attr(ans, "nHaplo") = nHaplo
  if (stabEval) {
    message("Computing stability")
    #   model fitting
    r2 = 1-norm(data-ans$haploStr%*%ans$haploFrq - ans$bias,"F")/norm(data-mean(data),"F")
    attr(ans, "R2") = r2
    #   weights stability
    omegaStab = vector('list', nBoot)
    for (j in 1:nBoot) {
      ind  = sample(1:n, round(n*fBoot))
      aux  = haploSep(data[ind,], nHaplo, stabEval=FALSE)
      indN = haploMatch(aux$haploStr, ans$haploStr[ind,])
      omegaStab[[j]] = aux$haploFrq[indN,]
    }
    lower = apply(simplify2array(omegaStab), 1:2,
                  function(x) quantile(x,0.025))
    upper = apply(simplify2array(omegaStab), 1:2,
                  function(x) quantile(x,0.975))
    attr(ans, "stabIntFrq") = data.frame(lowerBnd=lower, upperBnd=upper)
    attr(ans, "stabScoreFrq") = 1 - (upper - lower)
    #   structure stability
    fStab = vector('list', nBoot)
    for (j in 1:nBoot) {
      ind = sample(1:M, round(M*fBoot), replace=T)
      if (length(ind) < 2)
        stop(paste0(ncol(data), "*", fBoot, "cannot be smaller than 2"))
      fStab[[j]] = haploSep(data[,ind], nHaplo, stabEval=FALSE)$haploStr
    }
    resBoot = apply(simplify2array(fStab), 1:2,function(x) mean(x))
    scoreS = 1 - colMeans(abs(ans$haploStr - resBoot))
    attr(ans, "stabScoreStr") = scoreS
  }
  class(ans) <- c(class(ans), "haplo")
  return(ans)
}
