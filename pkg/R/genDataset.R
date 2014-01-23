#' @title genDataset
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description genDataset samples a random mixture with gaussian peak height.
#' @export
#' @usage genDataset(nC,popFreq)
#' @details Simple function for generating a random peak height mixture.
#' @param nC Number of contributors in model.
#' @param popFreq A list of allele frequencies for a given population.

genDataset = function(nC,popFreq,mu=10000,sd=1000) {
  require(gtools)
  nL<-length(popFreq)
  Mx=rdirichlet(1,rep(1,nC))  #true Mx for contributors
  refData <- list() 
  mixData <- list(adata=list(),hdata=list()) 
  for(i in 1:nL) {
   refData[[i]] <- list()
   mixA = numeric()
   mixH = numeric()
   for(s in 1:nC) {
    refData[[i]][[s]] = sample(names(popFreq[[i]]),size=2,prob=popFreq[[i]],replace=TRUE)
    mixA = c(mixA, refData[[i]][[s]] )
    mixH = c(mixH, abs(rnorm(2,mu*Mx[s],sd)))
   }
   agg=aggregate( mixH,by=list(mixA),sum)
   mixData$adata[[i]] = agg[,1]
   mixData$hdata[[i]] = agg[,2]
  }
  locs <- locnames(popFreq)
  names(mixData$adata) <- locs
  names(mixData$hdata) <- locs
  names(refData) <- locs
  return(list(trueMx=Mx,mixData=mixData,refData=refData))
}


