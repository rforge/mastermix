\name{calcDOdistr}
\alias{calcDOdistr}
\title{calcDOdistr}
\usage{
calcDOdistr(LRopt, sample, M = 1000, minS = 500, hdOnly = FALSE)
}
\arguments{
  \item{LRopt}{Model object returned by getLRopt function
  of GUI.}

  \item{sample}{is name of sample to be considered in the
  calculation.}

  \item{M}{is number of samples for each iteration.}

  \item{minS}{is number of minimum accepted samples.}

  \item{hdOnly}{is boolean wheter hp hypothesis is also
  considered.}
}
\value{
List with sample-elements for each hypothesis.
}
\description{
MCMC droppout sampler based on alleles in evidence.
}
\details{
calcDOdistr takes in models under Hp and Hd and samples
from the drop-out distribution based on total number of
alleles in evidence.
}
\author{
Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}

