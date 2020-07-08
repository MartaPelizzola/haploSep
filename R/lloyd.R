lloyd <- function (y, m, al, tsh, omega, relax, bias, weight,
                   lossFun, z, biasEst) {
  M <- ncol(y)
  n <- nrow(y)
  omdiff <- tsh
  if (missing(biasEst) || is.null(biasEst)) {
    biasV <- rep(0, M)
  } else {
    biasV <- biasEst
  }

  maxIt   <- 100 # 10e2
  countIt <- 0
  Al <- as.matrix(expand.grid(rep(list(al), m)))
  while (omdiff >= tsh && countIt <= maxIt) {
    val <- Al %*% omega[1:m,]
    if (lossFun == "LS") {
      pi <- sapply(1:nrow(y), function(x)
        which.min(colSums(weight[x,]*(y[x,] - biasV - t(val))^2)))
      ind <- sapply(1:m, function(x)
        which.min(colSums((t(y) - biasV - omega[x,])^2)))
      pi[ind] <- length(al)^(0:(m-1))+1
    } else {
      stop("It only supports least squares as loss function for now!")
    }
    
    if (!all(weight == 1) & bias){
      warning("The combination of bias =  True and weights is not available for now! The weight matrix will not be used here.")
      weight  = matrix(1,n,M)
    }
    omegaOld <- omega
    if (lossFun == "LS") {
      if (bias) {
        A = cbind(Al[pi,], rep(1, n))
        G = rbind( diag(c(rep(1,m),0)),-diag(c(rep(1,m), 0)))
        H = c(rep(0, m+1 ),c(rep(-1,m)), 0)
      } else {
        A = Al[pi,]
        G = rbind(diag(m),-diag(m))
        H = c(rep(0,m),rep(-1,m))
      }
      if (relax) {
        G = rbind(G, rep(-1,ncol(A)))
        H = c(H, -1)
        E = NULL
        F = NULL
      } else {
        E = rep(1,ncol(A))
        F = 1
      }
      if (all(weight == 1)) {
        aux <- sapply(1:M, function(x) lsei(A = A, B = y[,x], E = E, F = F,
                                            G = G, H = H, type = 2)$X)
        omega <- aux[1:m, ]
        if (bias)
          biasV <- aux[m+1, ]
      } else {
        aux <- sapply(1:M, function(x) lsei(A = A, B = y[,x], E = E, F = F,
                                            G = G, H = H, Wa = weight[,x])$X)
        omega <- aux[1:m, ]
        if(bias)
          biasV <- aux[m+1, ]
      }
    } else {
      stop("It only supports least squares as loss function for now!")
    }
    omdiff  <- max(abs(omega- omegaOld))
    countIt <- countIt + 1
  }

  if (omdiff >= tsh)
    warning(paste("Lloyds algorithm not converged with omega diff =",
                  round(omdiff, digits = 6), "and threshold =", tsh))

  ord   <- order(rowSums(omega),decreasing=TRUE)
  omega <- omega[ord,]
  pi    <- sapply(pi, function(x) which(colSums((Al[x,ord] - t(Al))^2) == 0))

  haploStr = Al[pi,]

  colnames(haploStr) <- sprintf("Haplo%d", 1:ncol(haploStr))
  rownames(omega)    <- sprintf("Haplo%d", 1:nrow(omega))

  list(haploFrq = omega, haploStr = haploStr, bias = biasV)
}


