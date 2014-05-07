###########
#Changelog#
###########
#18.04.13 - finish function and descriptions.
#09.10.13 -in helpfunction getPtab() g2 = as.numeric(unlist(strsplit(Q,''))) replaced with g2 = unlist(strsplit(Q,'')). 
#11.10.13 - bug: refData2[[k]]$adata[[k]][[i]] now given as refData2[[i]][[k]] 
#08.01.14 - We remove the "greedy"-method, uses OLS model for ind. loci-search and GLS for simultanously fit.
#23.01.14 - Changed from lociname to locinames
#24.01.14 - Changed order of input


#' @title deconvolve
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description deconvolve is a linear deconvolution procedure for STR DNA mixtures.
#' @export
#' @details The procedure optimizes the mixture proportion simultaneous with combined genotypes by assuming the STR response variation as normal distributed. The criterion for optimization is the error distance Mahalanobis Distance (MD) between the fitting model and observed responses.
#' 
#' Conditioning on referenced genotypes is possible. Selection of conditioned loci for each of the references may be specified. Unobserved alleles from references will be imputed as observed alleles with the input threshold as the quantitative information. Non-selected or empty observed loci will return NA as genotype combinations and not treated in model.
#' 
#' The search strategy is called keepElite which optimizes over all loci simultaniously by storing the 'eps' best fitted combinations during the search. The function also returns the optimized marginal result (each loci optimized).
#'
#' The covariance structures taking all loci into account assumes a compound symmetry structure which takes the number of alleles and peak heights into account (this ensures 'proportion of variance').
#' 
#' The user may choose whether combinations giving zero mixture propotion (gives overfitting model) for any contributors are accepted.
#' @param nC Number of contributors in model.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param refData Reference objects with list element [[i]][[s]] where list-element 's' is the reference index and the list-element 'i is the loci index where the qualitative data is stored as a length two vector.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
#' @param locsel_Mix Boolean-vector with Selected loci in mixData to deconvolve. locsel_Mix=NULL; accepts all loci.
#' @param locsel_Ref Boolean-matrix for specifying conditional loci (row) for each reference (column).locsel_Ref=NULL; accepts all loci.
#' @param eps Number of best combinations to keep during the search.
#' @param zeroMx boolean of allowing zero mixture proportion as an estimate for any contributors.
#' @param threshT Imputet quantitative value when conditioned reference alleles are non-observed.
#' @param verbose Boolean for whether in-process information should be printed
#' @return Optimized deconvolution model object.
#' \item{simpleList}{Table of loci independent optimizations. Uses independent covariance structure.}
#' \item{pList}{Resultlist of optimized combinations, mixture proportions and error-distances (MD).}
#' \item{locinames}{Name of loci in mixData}
#' \item{result1}{Tabled optimized results in format 1.}
#' \item{result2}{Tabled optimized results in format 2.}
#' \item{data}{All data used as input in analysis.}
#' \item{options}{Input parameters used in analysis.}
#' @references Tvedebrink,T, et.al.(2012). Identifying contributors of DNA mixtures by means of quantitative information of STR typing. Journal of Computational Biology, 19(7),887-902.
#' @keywords deconvolution, optimization

deconvolve = function(nC,mixData,refData=NULL,condOrder=NULL,locsel_Mix=NULL,locsel_Ref=NULL,eps=100,zeroMx=FALSE,threshT=50,verbose=FALSE) {
 require(MASS)
 require(gtools)
 #ERROR HANDLE:
 if(is.null(mixData$adata) | is.null(mixData$hdata)) { print('Missing mix-profile'); return(NULL) }
 Sa=length(mixData$adata)
 Sh=length(mixData$hdata)
 if(Sa!=Sh) { print('Wrong data format'); return(NULL) }
 nL = length(mixData$adata)
 if(!is.null(refData)) {
  if(length(refData)!=nL) { print("Wrong number of loci in refData"); return(NULL) }
  nR <- length(refData[[1]]) #number of refs
  if(any(sapply(refData,length)!=nR)) { print("Wrong number of refs for some loci"); return(NULL) }
  if(is.null(condOrder)) { print('CondOrder must be given');  return(NULL) }
  if(length(condOrder)!=nR) { print('Wrong length of condOrder given');return(NULL) }
  if(!is.null(locsel_Ref) && (nR!=ncol(locsel_Ref) || nL!=nrow(locsel_Ref))) { print('Wrong data format in locsel_Ref'); return(NULL) }
 }
 if(!is.null(locsel_Mix) && length(mixData$adata)!=length(locsel_Mix)) { print('Wrong data format in locsel_Mix'); return(NULL) }

 #insert lociname if not already given
 locinames = names(mixData$adata)
 if(is.null(locinames)) locinames = paste("Loci",1:length(mixData$adata),sep="")

 #Order-size too large
 if(!is.null(condOrder)) {
    msg = NULL
    if(is.null(refData)) msg <- paste('refData is missing')
    if(sum(condOrder>0)>nC) msg <- paste('You specified too many references',sep='')
    if(any(condOrder > nC)) msg <- paste('You specified too large number of condition-order',sep='')
    if(any(table(condOrder[condOrder>0])>1) ) msg <- paste('You specified atleast two references with same condition order.',sep='')
    if(!is.null(msg)) {
     stop(title='Conflict message', message=msg,icon='error')
     return(NULL) 
    }
 }
 
 #Modify data for deconvolution based on options. 
 mixData2 = mixData #data to modify
 refData2 = refData #data to modify
 #1) Unpresent loci to set threshT. This will increase number of alleles in data
 if(!is.null(condOrder)) {
  #goes through each cond. references (for each loci) to check for non-observed allele:
  for(i in 1:nL) {
   for(k in 1:length(condOrder)) {
    if(condOrder[k]==0 | length(refData2[[i]][[k]])==0) next
    aref = refData2[[i]][[k]]
    anew = aref[ !aref%in%mixData2$adata[[i]] ] #get new alleles not in data
    if(length(anew)>0) { 
     mixData2$adata[[i]] = c(mixData2$adata[[i]],anew) #add new alleles
     mixData2$hdata[[i]] = c(mixData2$hdata[[i]],rep(threshT,length(anew))) #add new alleles
     print(paste('WARNING: At locus ',locinames[i],', the allele(s) ',paste(anew,collapse="/",sep=""),' was added  with threshold height ',threshT,sep=''))
    }
   }
  }
 } 
 #2) Unchecked loci are set to numeric(). They will have length 0 and will be unselected.
 for(i in 1:nL) {
  if(!is.null(locsel_Ref)) {
   for(s in 1:ncol(locsel_Ref)) {
    if(!locsel_Ref[i,s]) { #check if loc i for reference j is checked
     refData2[[i]][[s]] = numeric()
    }
   }
  }  
  if(!is.null(locsel_Mix)) {
   if(!locsel_Mix[i]) { #loc.nr must be in selected loci
    mixData2$adata[[i]] = numeric()
    mixData2$hdata[[i]] = numeric()
   }
  }    
 }

 ######################### 
 #GENERIC HELP FUNCTIONS:#
 #########################
 #function for getting g-inverse covariance matrices
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

 fromAscii = function(x) { #takes vector with signs giving back integers 1 to 62
  asciRange = c(48:57,65:90,97:122) #max range of 62. This is type constant
  signs = sapply( as.raw(asciRange) , FUN = rawToChar)
  y = rep(NA,length(x))
  for(i in 1:length(y)) y[i] = grep(x[i],signs,fixed=TRUE,ignore.case = FALSE)
  return( y ) #distinguish small and large letters!
 }

 #function that get propotion of a vector
 getProps = function(h) {
  return(h/sum(h[!is.na(h)]))
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

 #Function returning allele names
 #numvec is vector with numbers
 getAnames = function(numvec,names) {
  names = names[!is.na(names)]
  nc = length(numvec) 
  ret = rep(NA,length(numvec))
  for(i in 1:nc) {
    ret[i] = paste( names[ fromAscii( unlist(strsplit(numvec[i],'')) ) ] ,collapse='/')
  }
  return(ret)
 }

 #function that renames the combinations in CC in order of Mx (for each combination) given (from major...)
 #returns list of contributors
 getCClist <-  function(CC,Alist,Mx,ordered=FALSE) {
  #Alist must be a list of locinames
  #CC is on form '12/22' '13/24' ...
  nL <- length(Alist)  
  nC <- length(unlist(strsplit(CC[1],'/')))
  #Need to init
  CClist <- list() #will have length #contributors
  for(k in 1:nC) { 
   CClist[[k]] <- list()
   for(i in 1:nL) {
    CClist[[k]][[i]] <- numeric(0)
   }
  }
  Mx_out = Mx
  for(k in 1:nC) {
   for(i in 1:nL) {
    for(l in 1:nrow(CC)) { #for each combination
     c=k #index of contributor
     if(ordered) {
      c <- order(Mx[l,],decreasing=TRUE)[k] #order by mixture prop
      Mx_out[l,] = sort(Mx[l,],decreasing=TRUE) #order Mx
     }
     subA <- fromAscii(unlist(strsplit(unlist(strsplit(CC[l,i],'/'))[k],'')))
     CClist[[c]][[i]] <- rbind(CClist[[c]][[i]],Alist[[i]][subA])
    }
   }
  }
  return(CClist=list(CC=CClist,Mx=Mx_out))
 }

 #Function filtering out the most probable combinations (used for MasterMix 1)
 #for a locus l in 1:L
 #ordering from most to less contribution
 getProbableList = function(l,Aname,MD,Mx,CC,eps,sortind,rel=FALSE) {
  #rel is mixuture proportion given relative to minor contributer
  #MD is Mahalanobis distance in normal- model
  #if eps=0, return the set with lowest MD
  #if eps<0, return constant number given
  #sortind is index to sort Mx of.
  if(all(is.na(MD))) return(NULL) #if all MD NA
  if(eps>=0) {
   indcrit = MD<=eps
   if(!any(indcrit,na.rm=TRUE)) { #if none choosen
    indcrit = min(MD,na.rm=TRUE) == MD
   }
  }
  if(eps<0) {
   nret = min(abs(eps),length(MD)-sum(is.na(MD))) #number to return for this loci. Difference with  umber of NA-values in MD
   indcrit = rep(FALSE,length(MD))
   indcrit[order(MD,decreasing=FALSE)[1:nret]]=TRUE #indices to use
  }
  indcrit[is.na(indcrit)] = FALSE
  if(!sum(indcrit)) return(NULL) #if no columns in CC, Return NULL - value
  tab = numeric(0)
  for(j in 1:sum(indcrit)) {
   if(sum(indcrit)==1) { cc = t(matrix(CC[indcrit,]))
  } else {  cc = CC[indcrit,][j,] }
   tab = rbind(tab,getAnames(cc,Aname))
  }
  #contribution values used to order arrangement of combinations
  mx = Mx[indcrit,]
  if(is.null(dim(mx))) mx = t(matrix(mx))
  ratios = mx
  for(j in 1:nrow(mx)) {
   if(rel) ratios[j,] = round(mx[j,]/min(mx[j,]),3)
   if(!rel) ratios[j,] = round(mx[j,],3)
   ordind = order(ratios[j,sortind],decreasing=TRUE) #select only those to sort
   ratios[j,sortind] = ratios[j,sortind[ordind]]
   tab[j,sortind] = tab[j,sortind[ordind]]
  }
  tab = cbind(l,tab,ratios,round(MD[indcrit],3))
  return(tab)
 }

 ########################## 
 #GENERIC MODEL FUNCTIONS:#
 ##########################

 #Function for generalized least squares estimatation (one locus)
 glinregfit <- function(Y,X,O,invW) {
  #Y is a datavector
  #X is covariance matrix
  #F is offset vector
  #invW is the inverse of correlation matrix
  XtWinv <- t(X)%*%invW
  mhat <- ginv(XtWinv%*%X)%*%XtWinv%*%(Y-O)
  MD <- t(Y-O)%*%invW%*%(Y-O) - 2*t(mhat)%*%XtWinv%*%(Y-O) + t(mhat)%*%XtWinv%*%X%*%mhat
#  Yhat <- X%*%mhat + O #predicted Y-respons from model.
  mhat <- c(mhat,1-sum(mhat)) #estimate of mhat
  return(list(mhat=mhat,MD = MD,Yhat=NULL))
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
 # Yhat <- list() #list for predicting responses
 # for(i in 1:length(Ylist)) { #for each elements in list:
 #  Yhat[[i]] <- Xlist[[i]]%*%mhat + Olist[[i]] #predicted response for loci i
 # }
  MD <- S_YtWinvY - 2*t(mhat)%*%S_XtWinvY + t(mhat)%*%S_XtWinvX%*%mhat
  mhat <- c(mhat,1-sum(mhat)) #estimate of mhat
  return(list(mhat=mhat,MD = MD,Yhat=NULL))
 }


 ###########################################################
 #**********************SIMPLE SEARCH *********************#
 ###########################################################
 #threating each loci independently
 #Function that gives proposed Mx and MD for given number of alleles 
 mastermix_simple = function(mixData2,nC,eps=0,refData2=NULL,condOrder=NULL,zeroMx=FALSE) {
  #mixData: Uses only replicate 1 and contains alleles and heights
  #nC: number assumed contributors at loci
  #refData: list with combinations from reference profiles (in alleles)
  #      -supports a subset of Loci (but both alleles must be in Evidence)
  #condOrder - vector for the specified references. specifies their order.
  #eps: 0 gives best combinations, >0 gives limit of MD, <0 gives max number of combinations
  if(is.null(eps)) stop("eps was not specified!")
  if(eps<0 && round(eps)!=eps) stop("eps must be a negative integer!")
  Alist = list() #listed allele info
  Ylist = list() #listed height info
  nL <- length(mixData2$adata)
  for(i in 1:nL ) { #for each loci
   Alist[[i]] = mixData2$adata[[i]]
   Ylist[[i]] = mixData2$hdata[[i]]
  }
  nA = unlist(lapply(Ylist,length))
  n = sum(nA) #total number of alleles
  invWi =  getCovmod(Ylist)$invWi  #get inverse covariance matricees

  #Mx and MD depends on comb of contributors which is different from loci
  #Store as list, with each list-elemnt as a given loci
  Mx = list()
  MD = list()
  CC = list() #contribute CC.
 
  for(i in 1:nL) { #for each loci
   if(nA[i]==0) {  #skip if no data at allele
    Mx[[i]] = NA
    MD[[i]] = NA
    CC[[i]] = NA
    next
   }
   #Method 1: each refs will be connected to a combination
   contr_combs <- getContrCombs(Alist=Alist[[i]],nC=nC,symmetry=FALSE,refs=refData2[[i]],condOrder=condOrder) #no symmetry!

   #Note: Aim is to assign some probaility of this set 
   #list-element of Mx-proposed and Residual sum squares
   nCombs = dim(contr_combs)[1]
   Mx[[i]] = matrix(NA,nrow=nCombs,ncol=nC)
   MD[[i]] = rep(NA,nCombs)
   CC[[i]] = contr_combs #stores ascii-signs here
   Sy <- sum(Ylist[[i]]) #get sum
 
   #Calc MD for each allele-comb
   for(j in 1:nCombs) { 
    Qi = contr_combs[j,] 
    Pi = matrix(getPtab(Qi),ncol=nC,nrow=nA[i]) #assign matrix: getPtab counts number of equal signs
    Pitilde <- matrix(Pi[,-nC],ncol=nC-1) #(n_i X (nC-1)) 
    Pic <- as.matrix(Pi[,nC])
    X = Sy/2*(Pitilde - Pic%*%t(matrix(1,nrow=nC-1))) #covariate matrix
    O = Sy/2*Pic #offset vector
    lmfit = glinregfit(Ylist[[i]],X,O,invW=invWi[[i]])
    if(any(is.na(lmfit$mhat)) | any(lmfit$mhat<0)) next  #only valid mx considered
    if(!zeroMx && any(lmfit$mhat==0)) next #skip when zero-contribution is not allowed
    Mx[[i]][j,] = lmfit$mhat
    MD[[i]][j] =  lmfit$MD
   }  #end each contributors
  } #end :i
  #create selection-list: is ranked
  #Note: Here the combinations are ranked from largest to lowest contributor
  #and hence the symmetry of the contributors is not necessary!
  simpleList = matrix(ncol=2*nC+2,nrow=0)
  for(i in 1:nL) {
   sortind = 1:nC
   if(!is.null(condOrder)) sortind = sortind[-condOrder]
   sublist = getProbableList(i,Alist[[i]],MD[[i]],Mx[[i]],CC[[i]],eps,sortind)
   if(is.null(sublist)) next
   sublist = sublist[order(sublist[,3],decreasing=TRUE),]
   simpleList = rbind(simpleList,sublist) 
  }
  if(nrow(simpleList)==0) return(NULL)
  rownames(simpleList) = 1:nrow(simpleList)
  simpleList = as.data.frame(simpleList,stringsAsFactors=FALSE)
  names(simpleList) = c('Loc',paste('Geno',1:nC,sep=''),paste('Mx',1:nC,sep=''),'MD')
  return(simpleList)
 } #END SIMPLE SEARCH


 ##############################################################
 #**********************KeepElite Search *********************#
 ##############################################################
 #The algorithm merges loci with decreasing order of number of alleles and minimum peak height and keep the best local subset of combinations (extension of greedy)
  #Finds all symmetric combinations but restricts on decreasing Mx for 'non-condition positions'.
 #Structure:
 #B - bunch-info
 #M - map-info of merged bunches
 #U - Unpealed combinations (survivors)
 #Procedure:
 #Loops until length(B)==1
 #Store combinations (with calc. D and Mx) with respect to some criterion.
 #Merge to new Bunches, store bunch-Map and loop through new combinations.
 #Store D and Mx when all bunches are calculated.
 mastermix_keepelite = function(mixData2,nC,eps=100,refData2=NULL, condOrder=NULL,zeroMx=FALSE,verbose=FALSE) {
  #mixData: Uses only replicate 1 and contains alleles and heights
  #nC: number assumed contributors at loci
  #refData[[i]][[s]]: list with combinations from reference profiles (in alleles)
  #      -supports a subset of Loci (but both alleles must be in Evidence)
  #condOrder - vector for the specified references. specifies their order.
  #eps: version1:proportion of extracted candidates, version2:max number of combinations
  if(is.null(eps)) stop("eps was not specified!")
  if(eps<=0 | round(eps)!=eps) stop("eps must be a positive integer!")

  cc = 0 #counter of number of loci
  Alist = list() #listed allele info
  Ylist = list() #listed height info
  locinames <- names(mixData2$adata)
  for(i in 1:length(locinames) ) { #for each loci
   cc = cc + 1
   Alist[[cc]] = mixData2$adata[[i]]
   Ylist[[cc]] = mixData2$hdata[[i]]
  }
  nL = cc #number of loci
  nA = unlist(lapply(Ylist,length)) #number of allele
  F= function(x) { if(length(x)==0) { return(0) } else { return(min(x)) } }
  minH = unlist( lapply(Ylist,F) )#minimum allele peak height
  n = sum(nA) #total number of alleles
  YiSum = unlist(lapply(Ylist,sum))
  Wi_list =  getCovmod(Ylist,OLS=FALSE)  #Uses GLS
 
  #1) Search first through locus with most information:
  #combinations: 
  done=FALSE
  seqcount <- 1 #counter of sequence
  B <- list() #list of current bunches
  C <- list() #list of survived combinations from each previous bunch-set
  nL2 = NULL
  #Outer while going through all combinations
  while(!done) {
   #initialize each loci to bunch if first seq:
   if(seqcount==1) { 
    cc = 0 #number of non-empty loci
    for(i in 1:nL) {
     if(nA[i]==0) next #should skip if 0 alleles
     cc = cc + 1
     B[[cc]] <- i
    }
    nbunchs <- nL2 <- cc 
   }
#   print(B)

   ##########################
   #fit model for each bunch#
   ##########################
   for(b in 1:nbunchs) { 
    if(seqcount==1) { 
     #Symmetry necessary: Using those with decreasing Mx
     contr_combs <- getContrCombs(Alist[[B[[b]]]],nC,symmetry=TRUE,refs=refData2[[B[[b]]]],condOrder=condOrder)#get combinations for loci
     C[[b]] = list(contr_combs)
    } else if(b>1) {
     next #only modelfit for bunch 1 
    }
    nLinB <- length(C[[b]]) #number of loci in bunch
    nCombs <- nrow(C[[b]][[1]]) #same number of rows in each list-element
    Mx <- matrix(NA,ncol=nC,nrow=nCombs)
    MD <- rep(NA,nCombs)
    for(j in 1:nCombs) { #for each combinations
     Y <- list()
     X <- list() 
     O <- list()
     iW <- list() 
     for(i in 1:nLinB) { #for each loci in Bunch
      Pi <- matrix(getPtab(C[[b]][[i]][j,]),ncol=nC,nrow=nA[B[[b]][i]])
      Pitilde = as.matrix(Pi[,-nC]) #(n_i X (nC-1)) 
      Pic = as.matrix(Pi[,nC])
      onetilde = as.matrix(rep(1,nC-1))
      X[[i]] = YiSum[B[[b]][[i]]]/2*(Pitilde - Pic%*%t(onetilde)) #expl. matrix
      O[[i]] = YiSum[B[[b]][[i]]]/2*Pic #offset vector 
      Y[[i]] = Ylist[[ B[[b]][[i]] ]] #add to list
      iW[[i]] = Wi_list$invWi[[ B[[b]][[i]] ]] #add to list
     }
     lmfit = glinregfit2(Y,X,O,iW) #GLS
     if(any(is.na(lmfit$mhat)) | any(lmfit$mhat<0)) next #only valid mx considered
     if(!zeroMx && any(lmfit$mhat==0)) next #skip when zero-contribution is not allowed
     Mx[j,] = lmfit$mhat
     MD[j] = lmfit$MD
     if(verbose && (j%%1000==0)) print(paste('comb',j,'finished for bunch',b))
    } #end each combination
    #all combinations for bunch b calculated:

    #If number of loci in bunch is all
    if(nLinB==nL2) { 
     done=TRUE
     break
    }

    #################
    ###KEEP ELITES###
    #################
    #peal off obvious errors:
    #CRITERION:
    #Pre-peal for first loci: Keep only those with correct Mx-order
    keep <- !is.na(MD) #1) keep those not NA
    if(seqcount==1) { #if first round.
     sortind <- 1:nC
     if(!is.null(condOrder)) sortind <- sortind[-condOrder[condOrder>0]] #don't sort referenced
     if(length(sortind)>0) {
      for(j in 1:length(MD)) { #for each fitted combination
       if(!keep[j]) next #skip those with negative Mx
       if(!all(order(Mx[j,sortind],decreasing=TRUE)==(1:length(sortind)))) keep[j] = FALSE #dont keep if not sorted Mx for non-references
      }
     }
    }
    if(sum(keep)==0) keep <- rep(TRUE,length(keep)) #keep all comb if all failed
    if(sum(keep)>eps) { #if still too many left
     keep2 <- MD[keep]<=sort(MD[keep])[eps] #keep the maxC best
     keep[keep] <- keep2
    }
    #keep only a specific combinations for bunch b:
    for(i in 1:nLinB) {
     C[[b]][[i]] <- matrix(C[[b]][[i]][keep,],nrow=sum(keep))
    }
   } #end for each bunch
   if(!done) { #continue to merge if not done
    #######################
    #MERGING and Updating:#  
    #######################
    #1) Get number of combinations in each bunch:
    nCombs = rep(NA,nbunchs)
    for(b in 1:nbunchs) {
     nCombs[b] <- nrow(C[[b]][[1]])
    }
    newB = list() #new Bunchlist
    M = list() #list for mapping merged bunches
    #We rank loci with number of alleles and minimum height, like the greedy algorithm
    #get rank of non-empty bunches
    nA2 = nA[nA>0] 
    minH2 = minH[nA>0] 
    if(seqcount==1) {  
     rank = order(nA2,decreasing=TRUE) #get ordered rank
     rankval = sort(nA2,decreasing=TRUE) #get ordered rank
     for(s in max(nA2):1) { #empty loci are not included in rank
      subrank = order(minH2[which(nA2==s)],decreasing=TRUE) #Rank after minH
      Si <- which(rankval==s) 
      rank[Si] = rank[Si][subrank] #change order of ranks
     }
#     print(rank)
    } else {
     rank = 1:nbunchs #the loci are already ranked
    }
    nB = nbunchs-1
    newNCombs = rep(NA,nB) #new number of combinations
    #update Bunchlist
    for(b in 1:nB) {
      if(b==1) ins = 1:2
      if(b>1) ins = b + 1
      newNCombs[b] <- prod(nCombs[rank[ins]])
      tmpB = numeric() #merging Bunches
      for(c in 1:length(ins)) {
       tmpB <- c(tmpB,B[[rank[ins[c]]]])
      }
      newB[[b]] = tmpB
      M[[b]] <- rank[ins]
    } 
    #update C (combination list)
    newC = list() 
    for(i in 1:length(newB)) { #for each element for current
     newC[[i]] = list()
     #for each prev bunc-element (by looking at bunchmap)
     for(j in 1:length(newB[[i]])) { 
      newC[[i]][[j]] = matrix(NA,nrow=newNCombs[i],ncol=nC) #stored combination matrix
     }
    }
  
    #for each combination: assign combinations as C[[bunch]][[lociind]] 
    for(i in 1:length(M)) { #for each element for current
     nloci = 0
     for(j in 1:length(M[[i]])) nloci <- nloci+length(B[[M[[i]][j]]])
     reccomb = 1 #determines how rapid combinations are changed
     for(j in length(M[[i]]):1) { #for each bunch in M, but from last to first
       for(k in length(B[[M[[i]][j]]]):1) { #for each locis in bunch (prev.) 
        tC <- C[[M[[i]][j]]][[k]] #each loci in previous Bunch
        Csz <- nrow(tC)
        antadd <- newNCombs[i]/Csz/reccomb
        reptC <- matrix(ncol=nC,nrow=reccomb*Csz)
        for(l in 1:nC) reptC[,l] <- c(t(replicate(reccomb,tC[,l])))
        toadd <- t(matrix(c(replicate(antadd,t(reptC))),nrow=nC))
        newC[[i]][[nloci]] <- toadd #add combination
        nloci <- nloci - 1 
       }
       reccomb <- reccomb*Csz #store sizes
     } #end for each bunch
    } #end Map-list 
    if(verbose) print(paste('#Combinations:',paste(newNCombs,collapse=' ')))
    nbunchs <- nB
    B <- newB
    C <- newC
    seqcount <- seqcount + 1
   }
  } #end repeat
  B <- B[[1]]
  C <- C[[1]]
  #sort and unlist Comb-list
  rank <- order(MD)
  Mx <- matrix(Mx[rank,],ncol=nC)
  MD <- MD[rank]  
  CC = matrix(ncol=nL,nrow=length(MD)) 
  for(i in 1:length(C)) {
   tC = matrix(C[[i]][rank,],ncol=nC)
   CC[,B[i]] <- apply(tC,1,function(x) { paste(x,collapse='/') })
  }
  CC[is.na(CC)] = paste(rep('NA',nC),sep='',collapse='/')
  #Order by major/mid/minor.. and get a list!
  pList <- list()
  CClist = getCClist(CC,Alist,Mx,ordered=FALSE)
  pList$CC <- CClist$CC
  pList$Mx <- CClist$Mx
  pList$MD <- MD #get distance
  return(pList)
 } #END KeepElite
 

####################################################################################################

 #CONTINUE THIS FUNCTION:

 #only assumes one replication:
 nA = unlist(lapply(mixData2$adata,length)) #number of alleles of selected loci
 #need to check if number of unique reference samples extends the model
 if(!is.null(refData2) && !is.null(condOrder)) {
  for(i in 1:nL) {
   refA = numeric()
   for(k in 1:length(condOrder)) {
    if(condOrder[k]==0) next
    refA = c(refA, refData2[[i]][[k]] )
   }
   #get number of references on locus i:
   nR = sum(condOrder>0)
   if(!is.null(locsel_Ref)) nR = sum( locsel_Ref[i,] &  condOrder>0) #if unselected markers in ref.
   refA = unique(refA) #unique refererence alleles (uncselected are not counted)
   Ai = mixData2$adata[[i]]
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
    return(NULL)
 }
 deconvlist <- list()
 deconvlist$pList <- mastermix_keepelite(mixData2,nC,eps,refData2,condOrder,zeroMx,verbose)
 deconvlist$simpleList <- mastermix_simple(mixData2,nC,eps=Inf,refData2, condOrder,zeroMx)

 #Helpfunction:
 #Function that returns all combinations given in returned deconv object
 getCombs <- function(deconvlist) {
  combrank = deconvlist$pList$CC
  nC <- length(combrank)
  nL <- length(combrank[[1]])
  nJ <- nrow(combrank[[1]][[1]]) #number of suggested combinations
  combs = matrix(nrow=nJ,ncol=nL)
  for(i in 1:nL ) { 
   for(j in 1:nJ ) {
    gencomb <- numeric(0)
    for(k in 1:nC) {
     geno <- paste(combrank[[k]][[i]][j,],collapse='/')
     gencomb <- c(gencomb,geno)
    }
    combs[j,i] <- paste(gencomb,collapse=',')
   }
  }
  ln <- deconvlist$locinames
  mxs <- rep("",nJ)
  for(j in 1:nJ) {
   mxs[j] <-  paste(round(deconvlist$pList$Mx[j,],3),collapse=",",sep="")
  }
  combs <- cbind(combs, mxs ,round(deconvlist$pList$MD,3) )
  rownames(combs) <- paste("#",1:nJ,sep="")
  colnames(combs) <- c(ln,'Mx','MD') #paste(c(t(replicate(2,ln))),c('_1','_2'),sep='')
  return(combs)
 }

 #Function same as getCombs2 but now listing with contributor:
 getCombs2 <- function(deconvlist) {
  combrank = deconvlist$pList$CC
  Mx <- deconvlist$pList$Mx
  MD <- deconvlist$pList$MD
  nC <- length(combrank)
  nL <- length(combrank[[1]])
  nJ <- nrow(combrank[[1]][[1]]) #number of suggested combinations
  combs = numeric(0) #matrix(nrow=nrow(combrank[[1]][[1]]*nC),ncol=length(combrank[[1]]))
  for(j in 1:nJ ) {
   for(k in 1:nC ) {
    gencomb <- numeric(0)
    for(i in 1:nL) { 
     gencomb <- c(gencomb,paste(combrank[[k]][[i]][j,],collapse='/'))
    }
    gencomb <- c(gencomb,round(Mx[j,k],2),round(MD[j],2))
    combs <- rbind(combs,gencomb)
   }
  }
  ln <- deconvlist$locinames
  rownames(combs) <- paste('#',c(t(replicate(nC,1:nJ))),':C',1:nC,sep='')
  ln[nchar(ln)==1] <- paste('0',ln[nchar(ln)==1],sep='')
  colnames(combs) <- c(ln,'Mx','MD') #paste(c(t(replicate(2,ln))),c('_1','_2'),sep='')
  return(combs)
 }

 #add extra result info:
 deconvlist$locinames <- names(mixData2$adata) #store locinames
 deconvlist$result1 = getCombs(deconvlist) #get results
 deconvlist$result2 = getCombs2(deconvlist)
 deconvlist$data <- list(mixData=mixData2,refData=refData2,locsel_Mix=locsel_Mix,locsel_Ref=locsel_Ref)
 deconvlist$options <- list(condOrder = condOrder,nC=nC,eps=eps)
 return(deconvlist)
} 
#***********************END deconvolution**************************#

