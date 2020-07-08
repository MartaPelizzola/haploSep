# Estimate haplotype weights with combinatorial algorithm
#
# @param y Numeric observation matrix with nrow(y) being the number of SNP locations and ncol(y) the number of time points.
# @param m Integer value which gives the number of essential haplotypes for which haplotype structure and weights should be estimated from the mixture.
# @param al Numeric vector which gives the possible values in the haplotype structure. By defauls al = c(0,1) corresponding to allel 0 and 1.
# @param bias a logical value indicating whether bias term (constant over SNP locations but different for various time points) should be added.
# @param round a logical value indicating whether pre-clustered values should be rounded to alphabet range.
#
# @return A matrix of dimension m x ncol(y) with estimated haplotype weights and if bias = TRUE with an attribute "bias" of estimated bias term.
# @references Behr, M., Munk, A., (2017), Identifiability for blind source separation of multiple finite alphabet linear mixtures, IEEE Transactions on Information Theory, 63(9):5506-5517
#


estOmega <- function(y, m, al = c(0,1), bias = TRUE, round = TRUE){
  n = nrow(y)
  M = ncol(y)

  if (length(al)^m > n)
    stop("Number of clusters is higher than number of SNPs. Use 'random' initialization in haploSep.")

  clusters   = hclust(stats::dist(y))
  clusterCut = cutree(clusters, length(al)^m)
  val = matrix(nrow = length(al)^m, ncol = M)
  for (i in 1:(length(al)^m)){
    ind.c = which(clusterCut == i)
    if (length(ind.c) > 1){
      val[i, ] = apply(y[ind.c,], 2, mean)
    } else {
      val[i, ] = y[ind.c,]
    }
  }
  val = val[order(rowSums(abs(val))),]
  if (bias) {
    biasEst = val[1, ]
    val = t(val) - val[1, ]
    val = t(val)
  }
  if (round) {
    val = pmin(val, max(al))
    val = pmax(val, min(al))
  }
  omega = matrix(nrow = m, ncol = M)
  omega[1,] = (val[2,]- al[1])/(al[2] - al[1])
  m.run = 1

  for (i in 2:m) {
    if (i == m) {
      alM = cbind(as.matrix(expand.grid(rep(list(al), m.run))), rep(0, length(al)^m.run))
    } else {
      alM = as.matrix(expand.grid(rep(list(al), m.run + 1)))
    }
    omega.run = rbind(omega[1:m.run,], rep(1, M) - colSums(omega, na.rm = T))
    val.run   = alM %*% omega.run

    ind = numeric(0)
    for (j in 1:dim(val.run)[1]) {
      val.order = order(colSums(abs(val.run[j,] - t(val))))
      ind.new = val.order[!sapply(val.order, function(x) is.element(x,ind))][1]
      ind = c(ind, ind.new)
    }
    omega[m.run + 1, ] = (val[-ind, ][1,] - al[1])/(al[2] - al[1])
    m.run = m.run + 1
  }

  if (bias) {
    attr(omega, "bias") <- biasEst
  }

  return(omega)
}
