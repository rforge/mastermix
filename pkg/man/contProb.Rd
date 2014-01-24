\name{contProb}
\alias{contProb}
\title{contProb}
\usage{
  contProb(nC, mixData, popFreq, refData = NULL,
    condOrder = NULL, locsel_Mix = NULL, locsel_Ref = NULL,
    M = NULL, threshT = 50)
}
\arguments{
  \item{nC}{Number of assumed contributors in model.}

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

  \item{locsel_Mix}{Boolean-vector with Selected loci in
  mixData to deconvolve. locsel_Mix=NULL; accepts all
  loci.}

  \item{locsel_Ref}{Boolean-matrix for specifying
  conditional loci (row) for each reference
  (column).locsel_Ref=NULL; accepts all loci.}

  \item{M}{Number of samples used in the importance
  sampling.}

  \item{threshT}{Imputet quantitative value when
  conditioned reference alleles are non-observed.}
}
\value{
  Either exact or an approximated evidence-probability with
  respect to defined hypothesis.
}
\description{
  contProb calculates the probability of the STR DNA
  mixture given a linear deconvolution model
}
\details{
  The procedure are doing either a exact probability
  calculation or an approximation by importance sample.
  Note that the probabilies are unnormalized with respect
  to the determinant of the covariance weights. Therefore
  only the ratio of two probabilities makes sense.

  The user may choose whether combinations giving zero
  mixture propotion (gives overfitting model) for any
  contributors are accepted.

  Function calls procedure in c++ by using the package
  Armadillo

  See Vignette for more details.
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Tvedebrink,T, et.al.(2012). Identifying contributors of
  DNA mixtures by means of quantitative information of STR
  typing. Journal of Computational Biology, 19(7),887-902.
}
\keyword{DNA}
\keyword{Importance}
\keyword{LR,}
\keyword{mixture}
\keyword{Quantitative}
\keyword{sample,}

