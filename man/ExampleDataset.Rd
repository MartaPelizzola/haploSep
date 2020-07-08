\name{Y}
\alias{Y}
\docType{data}
\title{
Allele frequencies from a simulated Evolve and Resequence experiment
}
\description{
Simulated allele frequencies with pool sequencing noise. An Evolve and 
Resequence experiment is simulated for 150 generations and the allele 
frequency data are available every 10 generations starting at time point 0.
}
\format{
  A matrix with nrow(Y)=500 observations on ncol(Y)=16 variables.
}
\details{
This simulated data are generated from a set of founder haplotypes from 
"Genetic redundancy fuels polygenic adaptation in Drosophila", 
Barghi et al. (2019). Additionally, haplotype frequencies are simualated 
to mimic an Evolve and Resequence experiment lasting 150 generations. 
Then allele frequencies are computed as the product of the haplotype
structre and frequency matrix. Pool sequencing noise, as a binomial 
sampling step, is added to the allele frequencies to obtain this data set.
}
\examples{
data(ExampleDataset)
str(Y) 
matplot(t(Y), type = "l", ylim = c(0,1), ylab = "Allele frequencies", xlab = "Generations", xaxt = "n")
axis(1, at = 1:16, labels=seq(0,150,10))
}
