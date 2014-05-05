\name{mastermixTK}
\alias{mastermixTK}
\title{mastermixTK}
\usage{
mastermixTK(envirfile=NULL)
}
\arguments{
  \item{envirfile}{A optional file to previous saved project of the program}
}
\description{
mastermixTK is a GUI for the function mastermix for
performing linear deconvolution procedure of STR DNA
mixtures.
}
\details{
The function is a graphical layer for the function
mastermix. See ?mastermix for more information.

The structure of profiles is:
Data[['samplename']]$adata[['locusname']] for allele
information and Data[['samplename']]$hdata[['locusname']]
for allele height information (if any)
}
\author{
Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}

\keyword{database search,}
\keyword{deconvolution,}
\keyword{Weight of evidence}

