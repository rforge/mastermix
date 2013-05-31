\name{deconvolve}
\alias{deconvolve}
\title{deconvolve}
\usage{
  deconvolve(mixData, nC, method = 'greedy', model =
    'general', eps = NULL, locsel_Mix = NULL, refData =
    NULL, locsel_Ref = NULL, condOrder = NULL, zeroMx =
    FALSE, threshT = 50)
}
\arguments{
  \item{mixData}{Evidence object with list elements
  adata[[i]] and hdata[[i]]. Each element has a loci-list
  with list-element 'i' storing qualitative data in 'adata'
  and quantitative data in 'hdata'.}

  \item{nC}{Number of contributors in model.}

  \item{method}{Selected search strategy: 'simple',
  'greedy' or 'peeloff'.}

  \item{model}{Selected covariance structure:
  'independent', 'weighted' or 'general'.}

  \item{eps}{Input parameter for search strategies:
  'independent': eps=0; keeps best combinations (lowest
  MD). eps>0; keeps combinations giving MD greater or equal
  eps. eps<0; keeps the abs(eps) best local combinations.
  'greedy': eps not used. 'peeloff': eps>0; number of best
  combinations to memorize thorough the search.}

  \item{locsel_Mix}{Boolean-vector with Selected loci in
  mixData to deconvolve. locsel_Mix=NULL; accepts all
  loci.}

  \item{refData}{Reference objects with list element
  [[s]]$adata[[i]]. The list element has reference-list
  with list-element 's' having a loci-list adata with
  list-element 'i storing qualitative data.}

  \item{locsel_Ref}{Boolean-matrix for specifying
  conditional loci (row) for each reference
  (column).locsel_Ref=NULL; accepts all loci.}

  \item{condOrder}{Specify conditioning references from
  refData (must be consistent order). For instance
  condOrder=(0,2,1,0) means that we restrict the model such
  that Ref2 and Ref3 are respectively conditioned as 2.
  contributor and 1. contributor in the model.}

  \item{zeroMx}{boolean of allowing zero mixture proportion
  as an estimate for any contributors.}

  \item{threshT}{Imputet quantitative value when
  conditioned reference alleles are non-observed.}
}
\value{
  Optimized deconvolution model object.
  \item{simpleList}{Table of loci independent
  optimizations} \item{pList}{Resultlist of optimized
  combinations, mixture proportions and error-distances
  (MD).} \item{res1}{Tabled optimized results in format 1.}
  \item{res2}{Tabled optimized results in format 2.}
  \item{data}{All data used as input in analysis.}
  \item{options}{Input parameters used in analysis.}
  \item{locinames}{Names of loci used in analysis (same as
  mixData$locinames if specified).}
}
\description{
  deconvolve is a linear deconvolution procedure for STR
  DNA mixtures.
}
\details{
  The procedure optimizes(sub) the mixture proportion
  simultaneous with combined genotypes by assuming the STR
  response variation as normal distributed. The criterion
  for optimization is the error distance Mahalanobis
  Distance (MD) between the fitting model and observed
  responses.

  Conditioning on referenced genotypes is possible.
  Selection of conditioned loci for each of the references
  may be specified. Unobserved alleles from references will
  be imputed as observed alleles with the input threshold
  as the quantitative information. Non-selected or empty
  observed loci will return NA as genotype combinations and
  not treated in model.

  The user may select between the search strategies
  'simple','greedy' and 'peeloff'. 'simple' optimizes each
  loci seperately, while the two latter optimizes over all
  loci simultaniously. The two latter differs in that
  thorough the search, 'greedy' memorizes the best
  combination, while the 'peeloff' memorizes multiple best
  combinations.

  The user may choose between different covariance
  structures; 'independent','weighted' or 'general' where
  the latter two satisfies the 'proportion of variance'.
  The latter also takes number of alleles into account and
  has the compound symmetry structure.

  The user may choose whether combinations giving zero
  mixture propotion (gives overfitting model) for any
  contributors are accepted.
}
\examples{
\dontrun{
genDataset = function(nC,nL,Aset = 1:6) {
 Mx=rdirichlet(1,rep(1,nC))  #true Mx for contributors
 refData <- list(adata=list())
 mixData <- list(adata=list())
 for(i in 1:nL) {
  mixA = numeric()
  mixH = numeric()
  for(s in 1:nC) {
   if(i==1) refData$adata[[s]] = list()
   refData$adata[[s]][[i]] = sample(Aset,size=2,replace=TRUE)
   mixA = c(mixA, refData$adata[[s]][[i]] )
   mixH = c(mixH, abs(rnorm(2,10000*Mx[s],50)))
  }
  agg=aggregate( mixH,by=list(mixA),sum)
  mixData$adata[[i]] = agg[,1]
  mixData$hdata[[i]] = agg[,2]
 }
 return(list(trueMx=Mx,mixData=mixData,refData=refData))
}
getRefs = function(refData) { #Function for returning reference genotype in a table
 R = numeric()
 nR = length(refData$adata)
for(k in 1:nR) {
 r = unlist(lapply( refData$adata[[k]], function(x) return(paste(x,collapse="/"))))
 R = cbind(R,r)
}
colnames(R) = paste("Ref",1:nR,sep="")
return(R)
}
###########
#EXAMPLE 1# Simple example with two contributors and 10 markers
###########
set.seed(0)
nC=2 #two contributors
nL=10 #10 markers
data = genDataset(nC,nL) #generate mixture data and reference profiles
mixData = data$mixData #observed mixture
deconv1 = deconvolve(mixData,nC,method='greedy',model='gen',zeroMx=FALSE)
deconv2 = deconvolve(mixData,nC,method='peeloff',model='gen',eps=50,zeroMx=FALSE)
print(deconv1$res2)
print(deconv2$res2[1:nC,])
print(t(rbind(getRefs(data$refData), signif(data$trueMx,3) )))
###########
#EXAMPLE 2# Advanced example with three contributors and 16 markers
###########
set.seed(0)
nC=3 #three contributors
nL=16 #16 markers
data = genDataset(nC,nL) #generate mixture data and reference profiles
mixData = data$mixData #observed mixture
refData = data$refData #true contributors
Rsel = c(3) #Selected References to condition on.
condOrder = rep(0,nC) #value zero means no condition of references
condOrder[Rsel]=1 #restrict Ref3 to position 1 in system.
refData$adata[[Rsel[1]]][[2]][1] = 7 #Ref 1 has unobserved in locus 2 (dropout i mix)
locsel_Mix = NULL #all loci considered in mixture
locsel_Ref = matrix(TRUE,ncol=length(condOrder),nrow=nL)
locsel_Ref[9,Rsel[1]] =FALSE #uncondition Ref 1 at locus 9 (ref from old kit)
deconv1 = deconvolve(mixData,nC,method="greedy",model="gen",eps=NULL,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
deconv2 = deconvolve(mixData,nC,method="peeloff",model="gen",eps=50,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
print(deconv1$res2)
print(deconv2$res2[1:nC,])
print(t(rbind(getRefs(refData), signif(data$trueMx,3) )))
###########
#EXAMPLE 3# Advanced example with four contributors and 10 markers
###########
set.seed(0)
nC=4 #four contributors
nL=10 #10 markers
data = genDataset(nC,nL) #generate mixture data and reference profiles
mixData = data$mixData #observed mixture
refData = data$refData #true contributors
Rsel = c(1,3) #Selected References to condition on (Ref1 and Ref3).
condOrder = rep(0,nC) #value zero means no condition of references
condOrder[Rsel]=c(1,2) #restrict Ref1 and Ref2 to position 1 and 2 in system.
locsel_Mix = NULL #no loci unselected
locsel_Ref = matrix(TRUE,ncol=length(condOrder),nrow=nL)
locsel_Ref[9,Rsel[1]] =FALSE #uncondition Ref 1 at locus 9
locsel_Ref[8,Rsel[2]] =FALSE #uncondition Ref 3 at locus 8
deconv1 = deconvolve(mixData,nC,method="greedy",model="gen",eps=NULL,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
deconv2 = deconvolve(mixData,nC,method="peeloff",model="gen",eps=50,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
print(deconv1$res2)
print(deconv2$res2[1:nC,])
print(t(rbind(getRefs(refData), signif(data$trueMx,3) )))
print(deconv2$res1[1:(50),])
}
}
\author{
  Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
}
\references{
  Tvedebrink,T, et.al.(2012). Identifying contributors of
  DNA mixtures by means of quantitative information of STR
  typing. Journal of Computational Biology, 19(7),887-902.
}
\keyword{deconvolution,}
\keyword{optimization}

