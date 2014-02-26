
#' @title contProbBayes
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contProbBayes calculates the probability of the STR DNA mixture given some assumed model by integrate out parameters.
#' @details The procedure are doing numerical integration to approximate the marginal probability by integrate over noisance parameters. Mixture proportions have flat prior.
#' 
#' Function calls procedure in c++ by using the package Armadillo
#'
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
#' @param COV Covariance structure of peak heights: "OLS", "WLS" or "GLS". Default is "WLS".
#' @param threshT The analytical threshold given. Used when considering allele drop-outs.
#' @param pTau Prior function for tau-parameter. Flat prior is default.
#' @param taumax Maximum range of tau-parameter. Default is 10000.
#' @param maxeval Maxumum number of evaluations in the interale function.
#' @return PE Aposterior probability of hypothesis (model) given evidence.
#' @references Tvedebrink,T, et.al.(2012). Identifying contributors of DNA mixtures by means of quantitative information of STR typing. Journal of Computational Biology, 19(7),887-902.
#' @keywords continous, Bayesian

#this model is introduced in Tvedebrink,T, et.al.(2012): But extended to bayesian model
contProbBayes = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,COV="WLS",threshT=50,pTau=function(x) { dgamma(1/x,1,0.001)}, taumax=10000, maxeval=10000){
 require(gtools)
 require(MASS) 
 require(cubature)

##########HELPFUNCTIONS############
 getCovmod = function(Ylist,COV="OLS") {
  invWi = list() #inverse is one to precalculate
  Wi=list()
  nA = unlist(lapply(Ylist,length))
  for(i in 1:length(nA)) {
   if(nA[i]==0) next
   if(COV=="OLS")  { Wi[[i]] <- diag(1,nA[i])
   } else if(COV=="WLS")  { Wi[[i]] <- diag(sum(Ylist[[i]]),nA[i])
   } else {
    Ci = diag(1,nA[i]) - matrix(1,nA[i],nA[i])/nA[i]
    Wi[[i]] = Ci%*%diag(Ylist[[i]],nA[i])%*%t(Ci) #pseudoinverse
   }
   invWi[[i]] <- ginv(Wi[[i]])
  }
  return(list(Wi=Wi,invWi=invWi))
 }
  #Function that get table of allele-counts from Q
  getPtab = function(Q) {
   g1 = c(t(replicate(2,1:length(Q))))
   g2 = unlist(strsplit(Q,''))
   gdat = rep(1,length(Q)*2)
   Ptab = tapply(X=gdat, INDEX=list(g1, g2), FUN=sum)
   Ptab[is.na(Ptab)] = 0 #must insert 0 manually
   return(t(Ptab))
  }
  #get population genotypes:
  getGlist <- function(popFreq) {
   locs <- names(popFreq)
   nL <- length(locs)
   Glist <- list()
   for(i in 1:nL) {
    G = t(as.matrix(expand.grid(rep(list(as.numeric(names(popFreq[[i]])),as.numeric(names(popFreq[[i]])) ))))) #one genotype per column
    keep = G[2,]>=G[1,] #unique genotypes 
    G <- G[,keep]  #store genotypes
    G <- matrix(as.character(G),nrow=2) #make string names again
    tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(popFreq[[i]]),as.numeric(popFreq[[i]]) )))))
    Gprob = exp(colSums(log(tmpP[,keep]))) #get allele probs
    ishet = G[1,]!=G[2,]
    Gprob[ishet] = 2*Gprob[ishet] #multiply with two to get heterozygote prob
    Glist[[locs[i]]] <- list(G=G,Gprob=Gprob)
   }
   return(Glist)
  }


  #function which checks if X==Y where rows and column are invariant. Use primenumbers
  isSameG = function(X,Y,rowInvariant=TRUE) {
   primtall = c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89, 101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211, 223,227,229,233)  
   if(any(dim(X)!=dim(Y))) return(FALSE) #wrong dimensions
   unX <- sort(unique(c(X)))
   unY <- sort(unique(c(Y)))
   if(length(unX)!=length(unY) || !all(unX==unY) ) return(FALSE) #different numbers
   C <- nrow(X)
   for(k in 1:length(unX)) { 
    X[X==unX[k]] <- primtall[k]
    Y[Y==unY[k]] <- primtall[k]
   }
   X <- matrix(as.integer(X),nrow=C,ncol=2)
   Y <- matrix(as.integer(Y),nrow=C,ncol=2)
   pX <- exp(rowSums(log(X)))
   pY <- exp(rowSums(log(Y)))
   if(!rowInvariant) return( all(pX==pY) )
   return(all( sort(pX)==sort(pY) ))
  }

  #function which returns combinations
  getQlist <- function(mixData,refData,condOrder) {
   Qcombs <- QcombsN <- list()
   nL <- length(mixData$adata)  
   for(i in 1:nL) {
    Alist <- mixData$adata[[i]] #==names(popFreq[[i]])?
    contr_combs <- Qcombs[[i]] <- getContrCombs(Alist,nC,symmetry=TRUE,refs=refData[[i]],condOrder,complete=TRUE) 
    ascii <- sort(unique(unlist(strsplit(contr_combs,"")))) 
    nCombs <- nrow(contr_combs)
    QcombsN[[i]] <- list()
    for(j in 1:nCombs) {  #for each genotype combination
     Qi = contr_combs[j,]
     cc <- lapply(strsplit(Qi,""),function(x) {
      tmp <- Alist[ascii%in%x]
      if(length(tmp)==1) tmp = rep(tmp,2)
      return(tmp)
     })
     QcombsN[[i]][[j]] <- t(matrix(unlist(cc),nrow=2))
     colnames(QcombsN[[i]][[j]]) <- paste("a",1:2,sep="")
     rownames(QcombsN[[i]][[j]]) <- paste("c",1:nC,sep="")
    }#end j
   } #end i
   return(list(Qcombs=Qcombs,QcombsN=QcombsN))
  } #end function

 #predefine all possible Q in Xlist,Olist. Also calculate genoProb
 getCombLists <- function(Qlist,Glist,nC,nAi,sY,condOrder) {
  maxQicomb = max(sapply(Qlist$Qcombs,nrow))
  pG <- matrix(nrow=nL,ncol=maxQicomb)
  Xlist <- Olist <- Plist <- matrix(list(),nrow=nL,ncol=maxQicomb)
  popG <- list()
  for(i in 1:nL) {
   for(j in 1:nrow(Qlist$Qcombs[[i]])) {  #for each genotype combination
      Qi = Qlist$Qcombs[[i]][j,]
      Tmp <- Qlist$QcombsN[[i]][[j]]
      pG[i,j] <- 1
      for(k in 1:nC) {
        if(!any(condOrder==k)) { #if no restriction
         a <- Glist[[i]]$G[1,]==Tmp[k,1] & Glist[[i]]$G[2,]==Tmp[k,2]
         b <- Glist[[i]]$G[1,]==Tmp[k,2] & Glist[[i]]$G[2,]==Tmp[k,1]
         pG[i,j] <- pG[i,j]*Glist[[i]]$Gprob[a | b] #genotype probabilities
        }
      }
      #factorial(length(Qi))/prod( factorial(table(Qi)) ) #permutation factor:
      Pi = matrix(getPtab(Qi),ncol=nC,nrow=nAi[i]) #assign matrix
      Pitilde <- matrix(Pi[,-nC],ncol=nC-1) #(n_i X (nC-1)) 
      Pic <- as.matrix(Pi[,nC]) 
      Xlist[[i,j]] <- sY[i]/2*(Pitilde - Pic%*%t(matrix(1,nrow=nC-1)))
      Olist[[i,j]] <- sY[i]/2*Pic #offset vector  
      Plist[[i,j]] <- Pi
   } #end comb:j
  } #end loci:i
  return(list(Xlist=Xlist,Olist=Olist,Plist=Plist,pG=pG))
 }

##########END HELPFUNCTIONS############

##########Start function here:############
  #Insert missing allele into response data:
  #taken from deconvolve()
   locinames <- names(popFreq)
   nL <- length(locinames)

   if (!is.null(condOrder)) {
    for (i in 1:nL) {
      for (k in 1:length(condOrder)) {
        if (condOrder[k] == 0 | length(refData[[i]][[k]]) == 0)   next
          aref = refData[[i]][[k]]
          anew = aref[!aref %in% mixData$adata[[i]]]
          if (length(anew) > 0) {
           mixData$adata[[i]] = c(mixData$adata[[i]],anew)
           mixData$hdata[[i]] = c(mixData$hdata[[i]],rep(threshT, length(anew)))
           print(paste("WARNING: At locus ",locinames[i],", the allele(s) ", paste(anew, collapse = "/", sep = ""), " was added  with threshold height ", threshT, sep = ""))
          }
        }
     }
  }
 nA = unlist(lapply(mixData$adata,length)) #number of alleles of selected loci
 #need to check if number of unique reference samples extends the model
 if(!is.null(refData) && !is.null(condOrder)) {
  nR = sum(condOrder>0)
  for(i in 1:nL) {
   refA = numeric()
   for(k in 1:length(condOrder)) {
    if(condOrder[k]==0) next
    refA = c(refA, refData[[i]][[k]] )
   }
   #get number of references on locus i:
   refA = unique(refA) #unique refererence alleles (uncselected are not counted)
   Ai = mixData$adata[[i]]
   leftA = Ai[!Ai%in%refA] #undescribed alleles for unknown
   if(length(leftA)>2*(nC-nR)) { 
    msg <- paste('For locus ',locinames[i],', number of unique allele left (after restriction) is ',length(leftA), ', while number of unknowns are ',nC-nR,'. Specify more unknowns or change reference conditioning.',sep='')
    stop(msg)
   }
  }
 }
 if(max(nA)>(2*nC)) {
    msg <- paste('Max alleles in a locus is ',max(nA),'. You must specify a greater number of contributors',sep='')
    stop(msg)
 }
 #END taken from deconvolve()

  Qlist <-  getQlist(mixData,refData,condOrder)
  Glist <- getGlist(popFreq)
  nAi <- sapply(mixData$adata,length)
  N <- sum(nAi)
  invWi = getCovmod(mixData$hdata,COV=COV)$invWi #used constant always
  Yvec <- unlist(mixData$hdata)
  sY = sapply(mixData$hdata,sum)
  lists <- getCombLists(Qlist,Glist,nC,nAi,sY,condOrder) 
  Ylist <- mixData$hdata
  Xlist <- lists$Xlist
  Olist <- lists$Olist
  nQi <- sapply(Qlist$Qcombs,nrow) #number of combs
  nQ <- sum(nQi)
#  print(paste("Number of combinations:",nQ,sep=""))
 
  cdfX <- cumsum( c(0,rep(nA*(nC-1),nQi)) )+1
  cdfO <- cumsum( c(0,rep(nA,nQi)) )+1 
  allX <- unlist(t(lists$Xlist)) #should be (C-1)x1 matrix
  allO <- unlist(t(lists$Olist))  #should be nx1 matrix
  allY <- unlist(Ylist)
  alliW <- unlist(invWi)
  pG <- c(t(lists$pG))
  pG <- pG[!is.na(pG)]
  PE = 1  

 prodlikY_omegaC <- function(theta) {   #call c++- function:
  val <- .C("calcLthetaC",as.numeric(PE),as.integer(nA), as.integer(nC), as.integer(nL), as.numeric(allY),  as.numeric(allX), as.numeric(allO),  as.numeric(alliW),  as.integer(cdfX), as.integer(cdfO),as.integer(nQi), as.numeric(pG), as.numeric(theta),PACKAGE="mastermix")[[1]]
  val <- val*pTau(theta[nC]); 
  return(val)
 }
 PEint <- adaptIntegrate(prodlikY_omegaC, lowerLimit = rep(0,nC), upperLimit = c(rep(1,nC-1),taumax),maxEval = maxeval )[[1]]
 return(PEint)
}

