
#' @title contProb
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contProb calculates the probability of the STR DNA mixture given a linear deconvolution model
#' @export
#' @details The procedure are doing either a exact probability calculation or an approximation by importance sample. Note that the probabilies are unnormalized with respect to the determinant of the covariance weights. Therefore only the ratio of two probabilities makes sense.
#'
#' The user may choose whether combinations giving zero mixture propotion (gives overfitting model) for any contributors are accepted.
#'
#' Function calls procedure in c++ by using the package Armadillo
#'
#' See Vignette for more details. 
#'
#' @param nC Number of assumed contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
#' @param locsel_Mix Boolean-vector with Selected loci in mixData to deconvolve. locsel_Mix=NULL; accepts all loci.
#' @param locsel_Ref Boolean-matrix for specifying conditional loci (row) for each reference (column).locsel_Ref=NULL; accepts all loci.
#' @param M Number of samples used in the importance sampling.
#' @param threshT Imputet quantitative value when conditioned reference alleles are non-observed.
#' @return Either exact or an approximated evidence-probability with respect to defined hypothesis.
#' @references Tvedebrink,T, et.al.(2012). Identifying contributors of DNA mixtures by means of quantitative information of STR typing. Journal of Computational Biology, 19(7),887-902.
#' @keywords Quantitative LR, Importance sample, DNA mixture

contProb = function(nC,mixData,popFreq,refData=NULL,condOrder=NULL,locsel_Mix=NULL,locsel_Ref=NULL,M=NULL,threshT=50){
 require(gtools)
 require(MASS)

##########HELPFUNCTIONS############
 getCovmod = function(Ylist,OLS=TRUE) {
  invWi = list() #inverse is one to precalculate
  Wi=list()
  nA = unlist(lapply(Ylist,length))
  for(i in 1:length(nA)) {
   if(nA[i]==0) next
   if(OLS)  { Wi[[i]] <- diag(1,nA[i])
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

#Function for generalized least squares estimatation (multiple loci)
 glinregfit2 <- function(Ylist,Xlist,Olist,invWlist) {
  #Ylist is a list of datavectors
  #Xlist is a list of covariance matrices
  #Flist is a list of offset vectors
  #invWlist is a list with inverse of correlation matrix CxC
 
  #precalcs of sums:
  S_YtWinvY <- 0 #1x1
  S_XtWinvY <- matrix(0,nrow=ncol(Xlist[[1]]),ncol=1) #(C-1)x1 matrix
  S_XtWinvX <- matrix(0,nrow=ncol(Xlist[[1]]),ncol=ncol(Xlist[[1]])) #(C-1)x(C-1) matrix
  for(i in 1:length(Ylist)) { #for each elements in list:
   Y <- Ylist[[i]]  
   X <- Xlist[[i]]  #varies with comb
   O <- Olist[[i]]  #varies with comb
   invW <- invWlist[[i]]
   XtWinv <- t(X)%*%invW
   S_YtWinvY <- S_YtWinvY + t(Y-O)%*%invW%*%(Y-O) #remember to differ from offset!
   S_XtWinvY <- S_XtWinvY + XtWinv%*%(Y-O)
   S_XtWinvX <- S_XtWinvX + XtWinv%*%X
  }  
  mhat <- ginv(S_XtWinvX)%*%S_XtWinvY
  MD <- S_YtWinvY - 2*t(mhat)%*%S_XtWinvY + t(mhat)%*%S_XtWinvX%*%mhat
  mhat <- c(mhat,1-sum(mhat)) #estimate of mhat
  return(list(mhat=mhat,MD = MD,Yhat=NULL))
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
  invWi = getCovmod(mixData$hdata,OLS=FALSE)$invWi #used constant always
  Yvec <- unlist(mixData$hdata)
  sY = sapply(mixData$hdata,sum)
  lists <- getCombLists(Qlist,Glist,nC,nAi,sY,condOrder) 
  Ylist <- mixData$hdata
  Xlist <- lists$Xlist
  Olist <- lists$Olist
  pG <- lists$pG
  nQiCombs <- sapply(Qlist$Qcombs,nrow) #number of combs

  #importance sampling:
  if(!is.null(M)) { 
  #deconvolve to obtain best combination:
   deconv <- deconvolve(nC=nC,mixData=mixData,refData=refData,condOrder=condOrder,locsel_Mix=NULL,locsel_Ref=NULL,eps=50,zeroMx=FALSE, threshT=threshT,verbose=FALSE)   #Obtain best combination:
   Qstar <- matrix(deconv$result2[1:nC,],nrow=nC) #get best combination
   Ghat <- list() #list for best genocombination
   for(i in 1:nL) Ghat[[i]] <- Qstar[,i]
  
   #1) get indices of best combination lists of X,O,..:
   BCGind <- rep(NA,nL)
   BXlist <- BOlist <- list()
   for(i in 1:nL) {
    tmpG = t(matrix(unlist(strsplit(Ghat[[i]],"/")),nrow=2))
    for(j in 1:length(Qlist$QcombsN[[i]])) {
     if(isSameG(tmpG,Qlist$QcombsN[[i]][[j]],rowInvariant=FALSE)) { #must have equal row-assignments
      BCGind[i] = j
      break; #found
     }
    }
    BXlist[[i]] <- Xlist[[i,BCGind[i]]]
    BOlist[[i]] <- Olist[[i,BCGind[i]]]
   }
   #2) calculate Qprob for different genocombinations
   LAprob <- function(MD,N,df=0) {
     sigmasq = MD/(N-df) #MLE of sigmasq
     return( (2*pi*sigmasq)^(-N/2)*exp(-0.5*(N-df)) ) #avoid |W| because it vanish in LR
   }
   LAvec <- matrix(nrow=nL,ncol=max(nQiCombs))
   for(i in 1:nL) {
    X <- BXlist #reset here 
    O <- BOlist #reset here 
    #go through all consistent combinations:
    for(j in 1:nQiCombs[i]) {  #for each genotype combination
       X[[i]] <- Xlist[[i,j]]
       O[[i]] <- Olist[[i,j]]
       lmfit = glinregfit2(Ylist,X,O,invWi) #GLS
       LAvec[i,j] <- LAprob(lmfit$MD,N,df=0)
    } #end comb:j
   } #end loci:i 
   qmargprob <- LAvec*pG  #multiply with genotype-probs here 
   for(i in 1:nL) {
    qmargprob[i,] <- qmargprob[i,]/sum(qmargprob[i,],na.rm=TRUE)
   }
   #3) Sampling from qmargprob:
   simX <- matrix(NA,ncol=nL,nrow=M)
   for(i in 1:nL) {
    isOK <- !is.na(qmargprob[i,])
    nG <- sum(!is.na(qmargprob[i,isOK]))
    simX[,i] <- sample(1:nG,M,replace=TRUE,prob=qmargprob[i,isOK])
   }
   simX = simX[do.call(order,as.data.frame(simX)),] # order the whole X-matrix (with elements>0):
  } #end if importance sampling

   #Prepare for C:
   Plist <- list()
   for(i in 1:nL) {
    Plist[[i]] <- numeric()
    for(j in 1:nQiCombs[i]) {
     Plist[[i]] <- rbind(Plist[[i]],lists$Plist[[i,j]]) #weights are Plist/Qmarg - values
    }
   }
   ULinvWi <- unlist(invWi)
   ULPlist <- unlist(Plist)
   if(!is.null(M)) pGvec <- c(t(lists$pG/qmargprob))
   if(is.null(M)) { 
    simX <- 0  
    M <- 0 #insert default
    pGvec <- c(t(lists$pG)) 
   }
   pGvec <- pGvec[!is.na(pGvec)]
   LA = 0
   retlist = .C("doAnalysisC",as.numeric(LA),as.integer(nAi),as.integer(nC),as.integer(nL),as.numeric(Yvec),as.numeric(ULinvWi),as.integer(nQiCombs),as.numeric(ULPlist),as.numeric(pGvec),as.integer(simX),as.integer(M),PACKAGE="mastermix")
   return(retlist[[1]])
} #end function


