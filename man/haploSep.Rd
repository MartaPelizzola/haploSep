% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/haploSep.R
\name{haploSep}
\alias{haploSep}
\title{Separates haplotypes}
\usage{
haploSep(
  data,
  nHaplo,
  stabEval = TRUE,
  bias = TRUE,
  weight = NULL,
  nBoot = 15,
  fBoot = 0.95
)
}
\arguments{
\item{data}{Numeric observation matrix with nrow(data) being the number of
SNP locations and ncol(data) the number of samples (e.g. time points).}

\item{nHaplo}{Integer value which gives the number of essential haplotypes
for which haplotype structure and frequency should be estimated from the
mixture. When missing, it is estimated via the SVD criterion.}

\item{stabEval}{logical. If \code{TRUE} (default), the stability analysis
of the reconstruction is carried out, otherwise only haplotype structure
and frequency are computed.}

\item{bias}{logical. It indicates whether bias term (constant over SNP
locations but different for various samples) should be added. By
default, it is \code{TRUE}}

\item{weight}{matrix of \code{ncol(weight)=ncol(data)} and
\code{nrow(weight)=nrow(data)} with weights for the observations least
squares optimization step. By default, equal weights are applied.}

\item{nBoot}{integer. It gives number of bootstrap repetitions for stability
score. By default, it takes value \code{15}. It is only needed when
\code{stabEval = TRUE}.}

\item{fBoot}{numeric. Its value lies between 0 and 1 which gives the
relative subsample size in the bootstrap runs. By default, it takes
value \code{0.95}. It is only needed when \code{stabEval = TRUE}.}
}
\value{
List with entries \code{haploFrq} and \code{haploStr},
which is an object from class haplo.

\itemize{
\item
\code{haploFrq} is a matrix with \code{nrow(haploFrq) = nHaplo}
and \code{ncol(haploFrq) = ncol(data)} which gives the estimated
 frequency of the estimated essential haplotypes.
\item
\code{haplotStr} is a matrix with \code{nrow(haploStr) = nrow(data)}
 and \code{ncol(haploStr) = nHaplo} and entries being either 0 or 1.
 It gives the estimated haplotype structure of essential haplotypes.
}

The returned list has an attribute "\code{nHaplo}", which is the number
of essential haplotypes.

If \code{stabEval = TRUE} the returned list has three more attributes
"\code{R2}", "\code{stabIntFrq}" and "\code{stabScoreStr}".

\itemize{
\item
The attribute "\code{R2}" is an inidicator of how good the model
fitting is, in a similar spirit as the \eqn{R^2} for linear models,
see \code{\link[stats]{lm}} and \code{\link[stats]{summary.lm}}.
\item
The attribute "\code{stabIntFrq}" provides a confidence envelope for
the estimated haplotype frequency. This is a data frame, containing
"\code{lowerBnd}" and "\code{upperBnd}" for each haplotype, which are
 0.025 and 0.975 quantiles for the bootstrap samples, respectively.
\item
The attribute "\code{stabScoreStr}" is the a numeric vector of length
\code{nHaplo} with values between 0 and 1.
}
}
\description{
Implementation of interative minimization algorithm for jointly estimating
structure of essential haplotypes as well as their relative proportion from
allele frequency matrix of mixture.
}
\examples{
 1. Reconstruct 5 haplotypes
 data(ExampleDataset)
 haploSep(data = Y, nHaplo = 5, stabEval = TRUE, bias = TRUE)
 
 Choose the number of haplotypes to be reconstructed with haploSelect
 data(ExampleDataset)
 m <- haploSelect(data = Y, bias = TRUE)
 haploSep(data = Y, nHaplo = m, stabEval = TRUE, bias = TRUE)
}
