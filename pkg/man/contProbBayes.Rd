\name{contProbBayes}
\alias{contProbBayes}
\title{contProbBayes}
\usage{
  contProbBayes(nC, mixData, popFreq, refData = NULL,
    condOrder = NULL, COV = "WLS", threshT = 50,
    pTau = function(x) { return(1) }, taumax=100 , maxeval = 10000)
}
\arguments{
  \item{nC}{Number of contributors in model.}

  \item{mixData}{Evidence object with list elements
  adata[[i]] and hdata[[i]]. Each element has a loci-list
  with list-element 'i' storing qualitative data in 'adata'
  and quantitative data in 'hdata'.}

  \item{popFreq}{A list of allele frequencies for a given
  population.}

  \item{refData}{Reference objects with list element
  [[s]]$adata[[i]]. The list element has reference-list
  with list-element 's' having a loci-list adata with
  list-element 'i storing qualitative data.}

  \item{condOrder}{Specify conditioning references from
  refData (must be consistent order). For instance
  condOrder=(0,2,1,0) means that we restrict the model such
  that Ref2 and Ref3 are respectively conditioned as 2.
  contributor and 1. contributor in the model.}

  \item{COV}{Covariance structure of peak heights: "OLS",
  "WLS" or "GLS". Default is "WLS".}

  \item{threshT}{The analytical threshold given. Used when
  considering allele drop-outs.}

  \item{pTau}{Prior function for tau-parameter. Flat prior
  is default.}

  \item{taumax}{Maximum range of tau-parameter. Default is
  10000.}

  \item{maxeval}{Maxumum number of evaluations in the
  interale function.}
}
\value{
  PE Aposterior probability of hypothesis (model) given
  evidence.
}
\description{
  contProbBayes calculates the probability of the STR DNA
  mixture given some assumed model by integrate out
  parameters.
}
\details{
  The procedure are doing numerical integration to
  approximate the marginal probability by integrate over
  noisance parameters. Mixture proportions have flat prior.

  Function calls procedure in c++ by using the package
  Armadillo
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Tvedebrink,T, et.al.(2012). Identifying contributors of
  DNA mixtures by means of quantitative information of STR
  typing. Journal of Computational Biology, 19(7),887-902.
}
\keyword{Bayesian}
\keyword{continous,}

