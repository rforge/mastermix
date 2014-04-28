\name{pvalueC}
\alias{pvalueC}
\title{pvalueC}
\usage{
pvalueC(LR.suspect, LR.table, P.table)
}
\arguments{
  \item{LR.suspect}{Observed LR to restrict sum on.}

  \item{LR.table}{matrix of all possible LR.}

  \item{P.table}{corresponding values to sum assigned in
  LR.table.}
}
\value{
Restriced sum from calculation.
}
\description{
Calculates Tippet-p-value in C.
}
\details{
This function calculates restricted sum for all values in
LR.table greater than LR.observed and summing corresponding
values in P.value.
}
\author{
Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}

