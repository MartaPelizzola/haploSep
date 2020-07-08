# Create random weight matrix
#
# Randomly create an m x M matrix with entries in (0,1) such that the M columns sum up to one.
#
# @param m Integer giving the number of rows.
# @param M Integer giving the number of columns.
#
# @return Matrix with m rows and M columns, entries in (0,1) such that the M columns sum up to one.
#

rOmega <- function(m,M){
  ans <- matrix(nrow=m, ncol=M)
  for (i in 1:M)
    ans[,i] <- (diff(sort(c(0,stats::runif(m-1),1))))
  ord <- order(rowSums(ans))
  ans <- ans[ord, ]
  as.matrix(ans, nrow = m, ncol = M)
}
