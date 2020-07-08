#' Match haplotype
#'
#' For a set of given haplotypes finds best matching haplotypes from a
#' canditate list with loss function being the number of non equal elements.
#'
#' @param haploStrA Matrix with given haplotype structure.
#' \code{nrow(haploStrA)} is number of SNP locations and \code{ncol(haploStrA)}
#' is number of estimated haploypes.
#' @param haploStrB Matrix of a candidate list of haplotypes.
#' \code{nrow(haploStrB)} is number of SNP locations and \code{ncol(haploStrB)}
#' is number of candidate haploypes.
#'
#' @return Integer vector with indexes of haploStrB which fit haploStrA best.
#'
#'
#'
#' @examples
#'  Compare the reconstruction with 3 dominant haplotypes to the one with 7 dominant haplotypes
#'  data(ExampleDataset)
#'  reconstruction1 <- haploSep(data = Y, nHaplo = 3, stabEval = FALSE, bias = TRUE)
#'  reconstruction2 <- haploSep(data = Y, nHaplo = 7, stabEval = FALSE, bias = TRUE)
#'  haploMatch(reconstruction1$haploStr, reconstruction2$haploStr)
#'  
#'  
#' @export


haploMatch <- function (haploStrA, haploStrB) {
  if (ncol(haploStrA) > ncol(haploStrB))
    stop('Input "haploStrA" should have no more columns than "haploStrB"!')
  if (nrow(haploStrA) != nrow(haploStrB))
    stop('Input "haploStrA" should have the same number of rows as "haploStrB"!')

  # via recursions
  # ndiff = -stats::cor(haploStrA,haploStrB)
  ndiff = t(haploStrA)%*%(1-haploStrB) + t(1-haploStrA)%*%haploStrB # exclusive or
  nh    = ncol(haploStrA)
  ord   = t(sapply(1:nh, function (x) order(ndiff[x,])[1:nh]))
  ndiff = t(sapply(1:nh, function (x) ndiff[x,ord[x,]]))
  idx   = rep(1,nh)
  dup   = duplicated(ord[1:nh+(idx-1)*nh])
  while (any(dup)) {
    fdl        = which.max(dup) # first duplicated location
    idx[1:fdl] = .rmOnlyDup(ndiff[1:fdl,], ord[1:fdl,], idx[1:fdl])
    dup        = duplicated(ord[1:nh+(idx-1)*nh])
  }
  id.oracle = ord[1:nh+(idx-1)*nh]
  return(id.oracle)
}

# Remove the only duplication
.rmOnlyDup <- function (ndiff, ord, idx, exList=NULL) {
  nh = length(idx)
  if (nh != nrow(ndiff) || nh != nrow(ord))
    stop('Dimension of input does not match!')
  if (sum(duplicated(ord[1:nh+(idx-1)*nh])) != 1)
    stop('There should be only one duplication!')
  if (nh == 2) {
    kopt = numeric(2)
    for (i in 1:2) {
      kopt[i] = 1
      while (ord[i,idx[i]+kopt[i]] %in% exList) {
        kopt[i] = kopt[i]+1
      }
    }
    li      = which.min(ndiff[1:2+(idx+kopt-1)*2])
    idx[li] = idx[li]+kopt[li]
    optIdx  = idx
  } else {
    dr   = which(duplicated(ord[1:nh+(idx-1)*nh])) # the right one in the duplicated pair
    err  = Inf
    i    = 0
    toGo = TRUE
    while (i < nh && toGo) {
      cIdx = idx
      cIdx[dr] = cIdx[dr]+i
      cExList  = exList
      while (ord[dr,cIdx[dr]] %in% cExList) {
        cIdx[dr] = cIdx[dr] + 1
        i = i + 1
      }
      cExList = c(cExList,ord[dr,cIdx[dr]])
      if (anyDuplicated(ord[1:nh+(cIdx-1)*nh])) {
        dl = which.max(ord[1:nh+(cIdx-1)*nh] == ord[dr,cIdx[dr]]) # the left one in the duplicated pair
        cIdx[dl] = cIdx[dl]+1
        while (ord[dl,cIdx[dl]] %in% cExList) {
          cIdx[dl] = cIdx[dl]+1
        }
        if (anyDuplicated(ord[1:nh+(cIdx-1)*nh])) {
          cIdx[-dr] = .rmOnlyDup(ndiff[-dr,],ord[-dr,],cIdx[-dr],cExList)
        }
      } else {
        toGo = FALSE
      }
      cErr = sum(ndiff[1:nh+(cIdx-1)*nh])
      if (cErr < err) {
        err    = cErr
        optIdx = cIdx
      }
      i = i+1
    }
  }
  optIdx
}


