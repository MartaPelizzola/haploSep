#' Simulate haplotype frequencies
#'
#' Given a matrix of haplotype structure simulates time series for haplotype 
#' frequencies under neutrality (drift only) or with selection. Returns 
#' haplotype frequency time series data, haplotype structure matrix, and the
#' matrix of allele frequency time series data computed from the haplotypes.
#' Sequencing  noise can be added to the allele frequency. In this case 
#' also a coverage matrix is returned.
#'
#' @param hp_str Haplotype structure matrix.  
#' \code{nrow(hp_str)} is number of SNP locations and \code{ncol(hp_str)}
#' is number of haploypes.
#' @param Ne One integer value corresponding to the effective population size.
#' @param tp A vecotr of integer values corresponding to the generations of
#' interest
#' @param meancov One integer value specifying the mean coverage to simulate
#' sequencing noise. If NULL the true allele frequencies are returned. 
#' @param hp0 Haplotype frequency matrix at generation 0. 
#' \code{length(hp0)} is number of haploypes. If NULL all haplotypes are 
#' assumed to start at equal frequency \code{1/ncol(hp_str)}
#' @param benef_all A vector of integer numbers corresponding to the rows
#' of hp_str with beneficial allele(s). If NULL neutral simulations are 
#' performed. 
#' @param s A vector of length \code{nrow(hp_str)} where 0 indicates that 
#' the corresponding allele is not beneficial, and selection coefficients 
#' are entered for each beneficial allele.
#' @param haploid Logical. TRUE for haploid simulations and FALSE for 
#' diploids.
#' 
#' @return List with allele frequency matrix, haplotype frequency matrix
#' and haplotype structure matrix. If \code{meancov!=NULL} the allele
#' frequency matrix contains noisy allele frequencies and an additional 
#' element corresponding to the coverage matrix is returned.
#'
#'
#'
#' @examples
#' hp_str <- matrix(sample(c(rep(0,20), rep(1,30)), 50), nrow = 10, ncol = 5)
#' Ne <- 300
#' tp <- seq(0,60,10)
#' meancov = 80
#' hp0 = c(0.1,0.3,0.2, 0.05, 0.35)
#' benef_all = c(1,3)
#' s = c(0.05,0,0.2,rep(0,7))
#' haploSimulate(hp_str, Ne, tp, meancov) #to simulate under neutrality. All haplotypes starting with equal frequency
#' haploSimulate(hp_str, Ne, tp, meancov, hp0) #to simulate under neutrality. Custom starting haplotype frequency.
#' haploSimulate(hp_str, Ne, tp, meancov,  benef_all, s) #to simulate with selection
#' 
#' @export




haploSimulate <- function(hp_str, Ne, tp, meancov=NULL, hp0 = NULL, benef_all=NULL, s=0, haploid=FALSE){
  if (missing(hp_str)){
    stop("The haplotype structure 'hp_str' is required to simulate haplotype frequencies")
  }
  if (missing(Ne)){
    stop("A value of the effective population size 'Ne' is required to simulate drift")
  }
  if (missing(tp)){
    stop("Time points of interest are required to simulate the haplotype frequencies")
  }
  if (!is.matrix(hp_str)){
    hp_str <- as.matrix(hp_str)
    warning("The haplotype structure was not of matrix type and it is now converted to a matrix")
  }
  if ((sum(hp_str == 0) + sum(hp_str==1))!=length(hp_str)){
    stop("All entries in the haplotype structure matrix 'hp_str' needs to be 0/1")
  }
  if (!(Ne==floor(Ne))){
    Ne <- as.integer(Ne)
    warning("Ne is transformed to the closest integer value")
  }
  if (!all(tp==floor(tp))){
    stop("Time points of interest need to be integer values")
  }
  if (is.null(hp0)){
    # Haplotype frequency constant at time point 0
    hp0 <- rep(1/ncol(hp_str),ncol(hp_str))
  }
  # Evolve the haplotypes
  max.len <- length(hp0)
  #forward in time simulations: initialize trajectory matrix
  hp_freq <- matrix(NA, ncol=length(tp), nrow=max.len, dimnames=list(c(), paste0("F", tp)))
  if(0 %in% tp){
    hp_freq[,"F0"] <- hp0
  }
  
  if(!haploid)
    Ne <- 2*Ne
  
  g <- 1
  p <- hp0
  wA <- 1+s 
  wa <- 1
  fitness_hp <- rep(1,max.len)
  hp_sel <- which(hp_str[benef_all,]==1)
  
  # simulate allele frequencies over time
  while(g <= max(tp)) {
    fitness_hp <- rep(1,max.len)
    if (sum(s!=0)>0){
      #apply selection
      for (i in 1:ncol(hp_str)){
        if(sum(hp_str[benef_all,i])>0){
          fitness_hp[i] <- fitness_hp[i]*sum(wA[benef_all][as.logical(hp_str[benef_all,i])])
        }
      }
      p <- p*fitness_hp
    }
    
    # apply drift
    if(!is.na(Ne))
      p <- rmultinom(1, Ne, p) / Ne
    
    if(g %in% tp){
      hp_freq[,paste0("F", g)] <- p
    }
    
    g <- g+1
  }  
  
  #True alelle frequency from the haplotypes
  afMat <- hp_str %*% hp_freq
  afMat <- pmin(pmax(afMat,0),1)
  if (!is.null(meancov)){
    covMat <- matrix(rpois(length(afMat), lambda=meancov), ncol=ncol(afMat))
    afMat_err <- rbinom(n = length(afMat), size = covMat, prob = afMat)/covMat #NOISY allele frequency
    return(list(afMat_err, hp_freq, hp_str, covMat))
  } else {
    return(list(afMat, hp_freq, hp_str))
  }
}
