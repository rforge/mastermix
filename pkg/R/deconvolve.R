###########
#Changelog#
###########
#18.04.13 - finish function and descriptions.

#' @title deconvolve
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description deconvolve is a linear deconvolution procedure for STR DNA mixtures.
#' @export
#' @usage deconvolve(mixData, nC, method = 'greedy', model = 'general', eps = NULL, locsel_Mix = NULL, refData = NULL, locsel_Ref = NULL, condOrder = NULL, zeroMx = FALSE, threshT = 50)
#' @details The procedure optimizes(sub) the mixture proportion simultaneous with combined genotypes by assuming the STR response variation as normal distributed. The criterion for optimization is the error distance Mahalanobis Distance (MD) between the fitting model and observed responses.
#' 
#' Conditioning on referenced genotypes is possible. Selection of conditioned loci for each of the references may be specified. Unobserved alleles from references will be imputed as observed alleles with the input threshold as the quantitative information. Non-selected or empty observed loci will return NA as genotype combinations and not treated in model.
#' 
#' The user may select between the search strategies 'simple','greedy' and 'peeloff'. 'simple' optimizes each loci seperately, while the two latter optimizes over all loci simultaniously. The two latter differs in that thorough the search, 'greedy' memorizes the best combination, while the 'peeloff' memorizes multiple best combinations. 
#'
#' The user may choose between different covariance structures; 'independent','weighted' or 'general' where the latter two satisfies the 'proportion of variance'. The latter also takes number of alleles into account and has the compound symmetry structure.
#' 
#' The user may choose whether combinations giving zero mixture propotion (gives overfitting model) for any contributors are accepted.
#' @param mixData Evidence object with list elements adata[[i]] and hdata[[i]]. Each element has a loci-list with list-element 'i' storing qualitative data in 'adata' and quantitative data in 'hdata'.
#' @param nC Number of contributors in model.
#' @param method Selected search strategy: 'simple', 'greedy' or 'peeloff'.
#' @param model Selected covariance structure: 'independent', 'weighted' or 'general'.
#' @param eps Input parameter for search strategies:  'independent': eps=0; keeps best combinations (lowest MD). eps>0; keeps combinations giving MD greater or equal eps. eps<0; keeps the abs(eps) best local combinations. 'greedy': eps not used. 'peeloff': eps>0; number of best combinations to memorize thorough the search.
#' @param locsel_Mix Boolean-vector with Selected loci in mixData to deconvolve. locsel_Mix=NULL; accepts all loci.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param locsel_Ref Boolean-matrix for specifying conditional loci (row) for each reference (column).locsel_Ref=NULL; accepts all loci.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model.
#' @param zeroMx boolean of allowing zero mixture proportion as an estimate for any contributors.
#' @param threshT Imputet quantitative value when conditioned reference alleles are non-observed.
#' @return Optimized deconvolution model object.
#' \item{simpleList}{Table of loci independent optimizations}
#' \item{pList}{Resultlist of optimized combinations, mixture proportions and error-distances (MD).}
#' \item{res1}{Tabled optimized results in format 1.}
#' \item{res2}{Tabled optimized results in format 2.}
#' \item{data}{All data used as input in analysis.}
#' \item{options}{Input parameters used in analysis.}
#' \item{locinames}{Names of loci used in analysis (same as mixData$locinames if specified).}
#' @references Tvedebrink,T, et.al.(2012). Identifying contributors of DNA mixtures by means of quantitative information of STR typing. Journal of Computational Biology, 19(7),887-902.
#' @keywords deconvolution, optimization
#' @examples
#' \dontrun{
#' genDataset = function(nC,nL,Aset = 1:6) {
#'  Mx=rdirichlet(1,rep(1,nC))  #true Mx for contributors
#'  refData <- list(adata=list()) 
#'  mixData <- list(adata=list()) 
#'  for(i in 1:nL) {
#'   mixA = numeric()
#'   mixH = numeric()
#'   for(s in 1:nC) {
#'    if(i==1) refData$adata[[s]] = list()
#'    refData$adata[[s]][[i]] = sample(Aset,size=2,replace=TRUE)
#'    mixA = c(mixA, refData$adata[[s]][[i]] )
#'    mixH = c(mixH, abs(rnorm(2,10000*Mx[s],50)))
#'   }
#'   agg=aggregate( mixH,by=list(mixA),sum)
#'   mixData$adata[[i]] = agg[,1]
#'   mixData$hdata[[i]] = agg[,2]
#'  }
#'  return(list(trueMx=Mx,mixData=mixData,refData=refData))
#' }
#' getRefs = function(refData) { #Function for returning reference genotype in a table
#'  R = numeric()
#'  nR = length(refData$adata)
#' for(k in 1:nR) {
#'  r = unlist(lapply( refData$adata[[k]], function(x) return(paste(x,collapse="/"))))
#'  R = cbind(R,r)
#' }
#' colnames(R) = paste("Ref",1:nR,sep="")
#' return(R)
#'}
#' ###########
#' #EXAMPLE 1# Simple example with two contributors and 10 markers
#' ###########
#' set.seed(0)
#' nC=2 #two contributors
#' nL=10 #10 markers
#' data = genDataset(nC,nL) #generate mixture data and reference profiles 
#' mixData = data$mixData #observed mixture
#' deconv1 = deconvolve(mixData,nC,method='greedy',model='gen',zeroMx=FALSE)
#' deconv2 = deconvolve(mixData,nC,method='peeloff',model='gen',eps=50,zeroMx=FALSE)
#' print(deconv1$res2)
#' print(deconv2$res2[1:nC,])
#' print(t(rbind(getRefs(data$refData), signif(data$trueMx,3) ))) 
#' ###########
#' #EXAMPLE 2# Advanced example with three contributors and 16 markers
#' ###########
#' set.seed(0)
#' nC=3 #three contributors
#' nL=16 #16 markers
#' data = genDataset(nC,nL) #generate mixture data and reference profiles 
#' mixData = data$mixData #observed mixture
#' refData = data$refData #true contributors
#' Rsel = c(3) #Selected References to condition on.
#' condOrder = rep(0,nC) #value zero means no condition of references
#' condOrder[Rsel]=1 #restrict Ref3 to position 1 in system.
#' refData$adata[[Rsel[1]]][[2]][1] = 7 #Ref 1 has unobserved in locus 2 (dropout i mix)
#' locsel_Mix = NULL #all loci considered in mixture 
#' locsel_Ref = matrix(TRUE,ncol=length(condOrder),nrow=nL)
#' locsel_Ref[9,Rsel[1]] =FALSE #uncondition Ref 1 at locus 9 (ref from old kit)
#' deconv1 = deconvolve(mixData,nC,method="greedy",model="gen",eps=NULL,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
#' deconv2 = deconvolve(mixData,nC,method="peeloff",model="gen",eps=50,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
#' print(deconv1$res2)
#' print(deconv2$res2[1:nC,])
#' print(t(rbind(getRefs(refData), signif(data$trueMx,3) ))) 
#' ###########
#' #EXAMPLE 3# Advanced example with four contributors and 10 markers
#' ###########
#' set.seed(0)
#' nC=4 #four contributors
#' nL=10 #10 markers
#' data = genDataset(nC,nL) #generate mixture data and reference profiles 
#' mixData = data$mixData #observed mixture
#' refData = data$refData #true contributors
#' Rsel = c(1,3) #Selected References to condition on (Ref1 and Ref3).
#' condOrder = rep(0,nC) #value zero means no condition of references
#' condOrder[Rsel]=c(1,2) #restrict Ref1 and Ref2 to position 1 and 2 in system.
#' locsel_Mix = NULL #no loci unselected
#' locsel_Ref = matrix(TRUE,ncol=length(condOrder),nrow=nL)
#' locsel_Ref[9,Rsel[1]] =FALSE #uncondition Ref 1 at locus 9 
#' locsel_Ref[8,Rsel[2]] =FALSE #uncondition Ref 3 at locus 8 
#' deconv1 = deconvolve(mixData,nC,method="greedy",model="gen",eps=NULL,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
#' deconv2 = deconvolve(mixData,nC,method="peeloff",model="gen",eps=50,locsel_Mix,refData,locsel_Ref,condOrder,zeroMx=FALSE)
#' print(deconv1$res2)
#' print(deconv2$res2[1:nC,])
#' print(t(rbind(getRefs(refData), signif(data$trueMx,3) ))) 
#' print(deconv2$res1[1:(50),])
#' }



deconvolve = function(mixData,nC,method='greedy',model='general',eps=NULL,locsel_Mix=NULL, refData=NULL,locsel_Ref=NULL,condOrder=NULL,zeroMx=FALSE,threshT=50) {
 require(MASS)
 require(gtools)
 #scaling of height is incorporated into model as Yi+.
 #ERROR HANDLE:

 if(is.null(mixData$adata) | is.null(mixData$hdata)) { print('Missing mix-profile'); return(NULL) }
 Sa=length(mixData$adata)
 Sh=length(mixData$hdata)
 if(Sa!=Sh) { print('Wrong data format'); return(NULL) }
 nL = length(mixData$adata)
 if(!is.null(refData$adata)) {
  nR = length(refData) #number of ref-profiles
  for(k in 1:nR) {
   if(length(refData[[k]]$adata)!=nL) { print('Wrong loci length in refData'); return(NULL) }
  }
 }
 if(!is.null(locsel_Mix) && length(mixData$adata)!=length(locsel_Mix)) { print('Wrong data format in locsel_Mix'); return(NULL) }
 if(!is.null(locsel_Ref) && length(refData)!=ncol(locsel_Ref)) { print('Wrong data format in locsel_Ref'); return(NULL) }
 if(!is.null(condOrder) && !is.null(refData) && length(refData)!=length(condOrder))  { print('Wrong length condOrder given'); return(NULL) }


 #insert lociname if not already given
 if(is.null(mixData$lociname)) mixData$lociname = names(mixData$adata)
 if(is.null(mixData$lociname)) mixData$lociname = paste("Loci",1:length(mixData$adata),sep="")

 if(any(length(grep(tolower(model),'ordinary')))) model=1
 if(any(length(grep(tolower(model),'weighted')))) model=2
 if(any(length(grep(tolower(model),'general')))) model=3


 #Order-size too large
 if(!is.null(condOrder)) {
    msg = NULL
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
 nL <- length(mixData2$lociname) #NUMBER OF LOCI
 #1) Unpresent loci to set threshT. This will increase number of alleles in data
 if(!is.null(condOrder)) {
  #goes through each cond. references to check for non-observed allele:
  for(k in 1:length(condOrder)) {
   for(i in 1:nL) {
    if(condOrder[k]==0 | length(refData2[[k]]$adata[[i]])==0) next
    aref = refData2[[k]]$adata[[i]]
    anew = aref[ !aref%in%mixData2$adata[[i]] ] #get new alleles not in data
    if(length(anew)>0) { 
     mixData2$adata[[i]] = c(mixData2$adata[[i]],anew) #add new alleles
     mixData2$hdata[[i]] = c(mixData2$hdata[[i]],rep(threshT,length(anew))) #add new alleles
     warning(paste('The allele(s) ',anew,' was added in locus ',i,' with threshold height ',threshT,sep=''))
    }
   }
  }
 } 
 #2) Unchecked loci are set to numeric(). They will have length 0 and will be unselected.
 for(i in 1:nL) {
  if(!is.null(locsel_Ref)) {
   for(s in 1:ncol(locsel_Ref)) {
    if(!locsel_Ref[i,s]) { #check if loc i for reference j is checked
     refData2[[s]]$adata[[i]] = numeric()
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
 getCovmod = function(Ylist,model) {
  invWi = list() #inverse is one to precalculate
  Wi=list()
  nA = unlist(lapply(Ylist,length))
  for(i in 1:length(nA)) {
   if(nA[i]==0) next
   if(model==1) Wi[[i]] <- diag(1,nA[i])
   if(model==2) Wi[[i]] <- diag(sum(Ylist[[i]]),nA[i])
   if(model==3) { 
    Ci = diag(1,nA[i]) - matrix(1,nA[i],nA[i])/nA[i]
    Wi[[i]] = Ci%*%diag(Ylist[[i]],nA[i])%*%t(Ci) #pseudoinverse
   }
   invWi[[i]] <- ginv(Wi[[i]])
  }
  return(list(Wi=Wi,invWi=invWi))
 }

 #Auxillary functions:
 toAscii = function(x) { #takes vector with integers 1 to 62
  asciRange = c(48:57,65:90,97:122) #max range of 62. This is type constant
  return( sapply( as.raw(asciRange[x]) , FUN = rawToChar) )
 }

 fromAscii = function(x) { #takes vector with signs giving back integers 1 to 62
  asciRange = c(48:57,65:90,97:122) #max range of 62. This is type constant
  signs = sapply( as.raw(asciRange) , FUN = rawToChar)
  y = rep(NA,length(x))
  for(i in 1:length(y)) y[i] = grep(x[i],signs,fixed=TRUE,ignore.case = FALSE)
  return( y ) #distinguish small and large letters!
 }
 
 #function that goes through perm-set of alleles and get genoypes using paste:
 permGenotype = function(permset) {
  permGen = numeric(0)
  for(i in 1:nrow(permset)) {
   permGen = c(permGen,paste(permset[i,],collapse=''))
  }
  return(permGen)
 }

 #function that get propotion of a vector
 getProps = function(h) {
  return(h/sum(h[!is.na(h)]))
 }
 
 #function taking in a set, number of samples and boolean as ordered
 getPerm = function(set,n,ordered) {
  if(!ordered) perm = combinations(length(set),n,v=set,repeats.allowed=TRUE)
  if(ordered) perm = permutations(length(set),n,v=set,repeats.allowed=TRUE)
  return(perm)
 }

 #method to get only allele with all contributors:
 removeCombsNotInA = function(combs,A) {
  temp = rep(TRUE,dim(combs)[1])
  bool = temp
  for(i in 1:length(temp)) temp[i] = paste(combs[i,],collapse='|')
  for(i in 1:length(A))   bool =   bool & grepl(A[i],temp)
  return(combs[bool,])
 }
 
 #Function that get table of allele-counts from Q
 getPtab = function(Q) {
  g1 = c(t(replicate(2,1:length(Q))))
  g2 = as.numeric(unlist(strsplit(Q,'')))
  gdat = rep(1,length(Q)*2)
  Ptab = tapply(X=gdat, INDEX=list(g1, g2), FUN=sum)
  Ptab[is.na(Ptab)] = 0 #must insert 0 manually
  return(t(Ptab))
 }

 #Function for getting allelecombinations for a given loci
 #When condition on reference, we need the order of contributors to be correct.
 #this is given by condOrder which indices contributors in a consistent order.
 #NB: If some ref-loci misses the loci is threated as unknown.
 getContr_combs <- function(Alist,i,nC,symmetry=FALSE,refs=NULL,condOrder=NULL) {
  #i is loci
  #symmetry=true gives all combinations
  #ref is list of alleles to condition on.
  #condOrder gives order of references in system. It must be same length as number of references in refs.
   #note: If a reference is unselected, it as numeric() in its loci
  Ai = toAscii(1:length(Alist[[i]]))
  singleCombs = getPerm(Ai,2,ordered=FALSE) #always two alleles pr. ind.
  if(any(Alist[[i]]%in%c('X','Y'))) {
   singleCombs = singleCombs[!rowSums(singleCombs==1)==2,]  #if amel: REMOVE YY-combination!
   if(is.null(dim(singleCombs))) singleCombs <- matrix(singleCombs,nrow=1)
  }
  genSet = permGenotype(singleCombs)
  #the alleles in ref must be a subset of alleles in Alist
  ncond = rep(FALSE,length(condOrder)) #boolean of vector wheter the condition is fine
  if(!is.null(refs)) { #if references given
   nRefs = length(refs) #number of references selected from import-frame
   condCC = rep(NA,nRefs) #vector of conditional combinations 
   for(k in 1:nRefs) { #for each reference we want to keep only combinations where these are included
    #check if all ref-alleles are in data A:
    #add the reference if both alleles are included:
    if(condOrder[k]==0) next #no condidtioning on the k-te reference
    condCC[k] = '' #init as no string
    if(length(refs[[k]]$adata[[i]])==0) next #continue if no reference there
    if(!all(refs[[k]]$adata[[i]]%in%Alist[[i]])) next #require both alleles to be in Alist
    CC = toAscii (which(Alist[[i]]%in%refs[[k]]$adata[[i]])) #get combinations 
    if(length(CC)==1) condCC[k] = paste(rep(CC,2),collapse='') #if hom
    if(length(CC)==2) condCC[k] = paste(sort(CC),collapse='') #if het, its sorted
    ncond[k] = TRUE #the conditional reference was accepted to use!
   }
  }
  multiCombs = getPerm(genSet,nC-sum(ncond),ordered=symmetry) #get unknown CC
  if(is.null(dim(multiCombs)))  return(t(as.matrix(multiCombs))) #return one combination
  finalcombs = matrix(NA,ncol=nC,nrow=nrow(multiCombs))
  #combine unknown CC with reference CC:
  #reference CC are ordered as in condOrder.
  if(!is.null(refs) & sum(ncond)>0) { #important to use the ncond here as check here!
   for(k in 1:length(ncond)) { #fill each condition
    if(!ncond[k]) next #skip ref if not included
    finalcombs[,condOrder[k]] = condCC[k] #or else insert conditional profile
   }
   finalcombs[,-condOrder[ncond]] = multiCombs #set in the randoms for those not condidtioned on
   multiCombs = finalcombs
  }
  contr_combs = removeCombsNotInA(multiCombs,Ai)
  if(is.null(dim(contr_combs)))  return(t(as.matrix(contr_combs)))
  return(contr_combs)
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

 getPriorityRec = function(i,tempvec,reclist,simpleList,nL) {
  #simpleList must be a table
  locind = which(i==simpleList$Loc) #get indices of loci
  if(length(locind)>0) {
   subComb =  simpleList[locind,]
   subComb =  subComb[order(subComb$MD),] 
   for(j in 1:dim(subComb)[1]) {
    tempvec = c(tempvec,locind[order(subComb$MD)][j])
    if(i==nL) { #if last loci
     reclist = rbind(reclist, tempvec) #store comb in the end here
    } else {
     reclist = getPriorityRec(i+1,tempvec,reclist,simpleList,nL) #receive returned list
    }
    tempvec = tempvec[-i] #remove i'te element when returned
    } #end loop
   } else {
    tempvec = c(tempvec,0)
   if(i==nL) { #if last loci
    reclist = rbind(reclist, tempvec) #store comb in the end here
   } else {
    reclist = getPriorityRec(i+1,tempvec,reclist,simpleList,nL) #receive returned list
   }
  }
  return(reclist) #return reclist here
 }

 #function returning pList for combinations in simpleList:
 #averaging Mx for multiple loci.
 getpList <- function(simpleList,Alist) {
  nL <- length(Alist)
  recList = getPriorityRec(i=1,tempvec=numeric(0),reclist=numeric(0),simpleList,nL)
  cind = grep('Geno',names(simpleList)) #columns of contrs
  sind = grep('MD',names(simpleList)) #column of MD
  mind = grep('Mx',names(simpleList)) #column of Mx
  nC <- length(mind)
  pList <- list()
  pList$CC <- list()
  pList$MD = rep(0,nrow(recList))
  pList$Mx = matrix(0,nrow(recList),ncol=nC)

  #NB: Mx is already sorted from greatest to lowest contributors
  #Mx is averaged over non-empty loci
  for(k in 1:nC) { #each contributor. 
   pList$CC[[k]] <- list()
   for(i in 1:nL) { #each loci
    a <- simpleList[recList[,i],cind[k]] #take out allele names
    if(length(a)>0) {
     pList$CC[[k]][[i]] <- t(matrix(unlist(strsplit(a,'/')),nrow=2))
    } else if(length(Alist[[i]])==1) {
     pList$CC[[k]][[i]] <- matrix(Alist[[i]],ncol=2,nrow=nrow(recList))
    } else {
     pList$CC[[k]][[i]] <- matrix(NA,ncol=2,nrow=nrow(recList))
    }
    if(k==1) {
     if(all(recList[,i]==0)) next
     pList$MD <- pList$MD + as.numeric(simpleList[recList[,i],sind])
     pList$Mx <- pList$Mx + matrix(as.numeric(as.matrix(simpleList[recList[,i],mind])),ncol=nC)
    }
   }
  }
  pList$Mx <- pList$Mx/sum(recList[1,]>0) 
  return(pList)
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
 #mastermix_simple: MastermMix1v2.R
 #threating each loci independently
 #Function that gives proposed Mx and MD for given number of alleles 
 mastermix_simple = function(mixData2,nC,eps=0,model=1,refData2=NULL, condOrder=NULL, pList=TRUE,zeroMx=FALSE) {
  #####MasterMix1-helpfunctions#######
  #Help function for recursive looking:
  #if nA=1, we know that it is
  #model={1,2,3} - OLS, WLS or GLS
  #mixData: Uses only replicate 1 and contains alleles and heights
  #nC: number assumed contributors at loci
  #refData: list with combinations from reference profiles (in alleles)
  #      -supports a subset of Loci (but both alleles must be in Evidence)
  #condOrder - vector for the specified references. specifies their order.
  #pList - boolean whether to carry out final list
  if(is.null(eps)) stop("eps was not specified!")
  if(eps<0 && round(eps)!=eps) stop("eps must be a negative integer!")

  cc = 0 #counter of number of loci 
  Alist = list() #listed allele info
  Ylist = list() #listed height info
  for(i in 1:length(mixData2$lociname) ) { #for each loci
   cc = cc + 1
   Alist[[cc]] = mixData2$adata[[i]]
   Ylist[[cc]] = mixData2$hdata[[i]]
  }
  nL = cc #number of loci
  nA = unlist(lapply(Ylist,length))
  n = sum(nA) #total number of alleles
  Wi_list =  getCovmod(Ylist,model)  #get inverse covariance matricees

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
   #No symmetry. Fix conditinional contributors
   contr_combs <- getContr_combs(Alist,i,nC,symmetry=FALSE,refData2,condOrder) #no symmetry!

   #Note: Aim is to assign some probaility of this set 
   #list-element of Mx-proposed and Residual sum squares
   nCombs = dim(contr_combs)[1]
   Mx[[i]] = matrix(NA,nrow=nCombs,ncol=nC)
   MD[[i]] = rep(NA,nCombs)
   CC[[i]] = contr_combs #stores ascii-signs here
   yi = Ylist[[i]]
 
   #Calc MD for each allele-comb
   for(j in 1:nCombs) { 
    Qi = contr_combs[j,] 
    Pi = matrix(getPtab(Qi),ncol=nC,nrow=nA[i]) #assign matrix: getPtab counts number of equal signs
    Pitilde <- matrix(Pi[,-nC],ncol=nC-1) #(n_i X (nC-1)) 
    Pic <- as.matrix(Pi[,nC])
    X = sum(yi)/2*(Pitilde - Pic%*%t(matrix(1,nrow=nC-1))) #covariate matrix
    O = sum(yi)/2*Pic #offset vector
    lmfit = glinregfit(yi,X,O,invW= Wi_list$invWi[[i]])
    if(any(is.na(lmfit$mhat)) | any(lmfit$mhat<0)) next  #only valid mx considered
    if(!zeroMx && any(lmfit$mhat==0)) next #skip when zero-contribution is not allowed
  
    #stored comb. requirement: Mx decreases + condOrder on right indices
    #if(!all(sort(lmfit$mhat,decreasing=TRUE)==lmfit$mhat)) next #
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
  #receive combinations and order them by ascending MD
  if(pList) {
   print('Receiving prioritylist...')
   pList <- getpList(simpleList,Alist) #get listed structure on contributors
   print('...DONE')
  }
  #get simple analysis in addition
  if(!is.infinite(eps)) simpleList = mastermix_simple(mixData2,nC,eps=Inf,model,refData2, condOrder, pList=FALSE,zeroMx)$simpleList #update simpleList if not inf given
  return(list(simpleList=simpleList,pList=pList))
 } #END SIMPLE SEARCH




 ###########################################################
 #**********************Greedy Search *********************#
 ###########################################################
 #The algorithm merges loci with decreasing order of number of alleles and minimum peak height and keep the local best combination
 #Note: Finds all symmetric combinations but restricts on decreasing Mx for 'non-condition positions'
 mastermix_greedy = function(mixData2,nC,eps=NULL,model=1,refData2=NULL, condOrder=NULL,zeroMx=FALSE) {
  cc = 0 #counter of number of loci
  Alist = list() #listed allele info
  Ylist = list() #listed height info
  for(i in 1:length(mixData2$lociname) ) { #for each loci
   cc = cc + 1
   Alist[[cc]] = mixData2$adata[[i]]
   Ylist[[cc]] = mixData2$hdata[[i]]
  }
  nL = cc #number of loci
  nA = unlist(lapply(Ylist,length)) #number of allele
  F= function(x) { if(length(x)==0) { return(0) } else { return(min(x)) } }
  minH = unlist( lapply(Ylist,F) )#minimum allele height
  n = sum(nA) #total number of alleles
  
  Wi_list =  getCovmod(Ylist,model)  #get inverse covariance matricees
  simpleList = list() #List for storing 
  YiSum = unlist(lapply(Ylist,sum))
  D <- Inf
  mhat <- rep(0,nC)
  loc <- rep(NA,nL)
 
  #Greedy algoritm:
  filled <- FALSE #visited all loci?
  Y <- list()
  W <- list()
  X <- list()
  O <- list()
  CC <- list()
  nonSymInd = list() #list of index in contr_combs to use. Calc. in first round

  repeat{ #repeat until tau is unchanged
   #values to accumulate
   anyBetter <- FALSE #boolean if any loci gave new updated value
   count <- 1 #loci visited-counter
   for(i in (2*nC):0) {  #note: goes down to zero, but zero are not calculated!
    Si <- which(nA==i)
    newrank <- order(minH[Si],decreasing=TRUE)  #order Si after minH: minumum height in loci
    Si <- Si[newrank] #update new rank
    nS <- length(Si)
    if(nS==0) next #if no observation for block: Go to next block!
    for(j in 1:nS) { #for each loci withing block
     loc[count] <- Si[j]
     if(nA[ Si[j]]==0) {
      CC[[count]]=  cbind(paste(rep('NA',nC),sep='',collapse='/')) #insert NA/NA in combination if no data
      count <- count + 1
      next #if no data in locus we skip
     }
     Y[[count]] = Ylist[[Si[j]]] #necessary to keep Y in own vector
    
     #Symmetry is used to get combinations which gives decreasing order of Mx
     #this is a necessary in the algorithm that contr. stays with relative same ranks of mx in OLS-model
     #order contributors by the order in refs 
     contr_combs <- getContr_combs(Alist,Si[j],nC,symmetry=TRUE,refData2,condOrder)#get combinations for loci. condOrder gives column of conditional references 
     #At first round we must find out which combinations are major/middle/minor
     #this is necessary to do since the algorithm is nesting
     #store what combinations to use further in ALGORITHM
     if(!filled) { #if first round
      nonSymInd[[count]] <- numeric()
      for(k in 1:nrow(contr_combs)) { 
       Qi = contr_combs[k,]
       Pi = matrix(getPtab(Qi),ncol=nC,nrow=nA[Si[j]]) #assign matrix
       Pitilde <- matrix(Pi[,-nC],ncol=nC-1) #(n_i X (nC-1)) 
       Pic <- as.matrix(Pi[,nC])
       X2 <- YiSum[Si[j]]/2*(Pitilde - Pic%*%t(matrix(1,nrow=nC-1)))
       O2 <- YiSum[Si[j]]/2*Pic #offset vector 
       lmfit = glinregfit(Y[[count]],X2,O2,invW=Wi_list$invWi[[Si[j]]])
       #Keep only valid estimates and decreasing proportion estimates: k-te combination is used for later
       if(any(is.na(lmfit$mhat)) | any(lmfit$mhat<0)) next #only valid mx considered
       if(!zeroMx && any(lmfit$mhat==0)) next #skip when zero-contribution is not allowed

       #if condOrder have order in in contr_combs, the first length(refs$adata) is fixed and the rest have decreasing Mx
       sortind <- 1:nC
       if(!is.null(condOrder)) sortind <- sortind[-which(condOrder>0)] #don't sort referenced
       #unconditional combinations are sorted by mx-hat. If they are decreasing, the combinations are stored.
       if(length(sortind)==0) nonSymInd[[count]] <- c(nonSymInd[[count]],k) #all are conditioned if sortind has length zero and are accpected.
       if(length(sortind)>0 && all(sort(lmfit$mhat[sortind],decreasing=TRUE)==lmfit$mhat[sortind])) nonSymInd[[count]] <- c(nonSymInd[[count]],k) #the combination was valid
      } #end k:
     } #end if not filled

     #Update cases with only when valid mhat is found
     if(length(nonSymInd[[count]])>0) contr_combs <- contr_combs[nonSymInd[[count]],] #update which combinations to use
     if(is.null(dim(contr_combs))) contr_combs = t(as.matrix(contr_combs))
     nCombs <- nrow(contr_combs)
     tempMx <- matrix(nrow=nCombs,ncol=nC)
     tempD <- rep(NA,nCombs)
     tempCC <- rep(NA,nCombs)
     tempX <- list()
     tempO <- list()
     tempW <- list()
 
     for(k in 1:nCombs) {  #for each genotype combination
      Qi = contr_combs[k,]
      Pi = matrix(getPtab(Qi),ncol=nC,nrow=nA[Si[j]]) #assign matrix
      Pitilde <- matrix(Pi[,-nC],ncol=nC-1) #(n_i X (nC-1)) 
      Pic <- as.matrix(Pi[,nC]) 
      X[[count]] <- tempX[[k]] <- YiSum[Si[j]]/2*(Pitilde - Pic%*%t(matrix(1,nrow=nC-1)))
      O[[count]] <-  tempO[[k]] <- YiSum[Si[j]]/2*Pic #offset vector 
      CC[[count]] <- tempCC[[k]] <- paste(Qi,collapse='/')
 
      #fit model
      W[[count]] <- tempW[[k]] <- Wi_list$invWi[[Si[j]]] #GLS
      lmfit = glinregfit2(Y,X,O,W) #GLS
      if(any(is.na(lmfit$mhat)) | any(lmfit$mhat<0)) next
      tempMx[k,] <- lmfit$mhat #store Mx-fit of all combinations
      tempD[k] <- lmfit$MD #store distance of all combinations
     } #end comb:k

     ch <- which.min(tempD)
     #remember to update storage-list with the local best!
     X[[count]] <- tempX[[ch]]
     O[[count]] <-  tempO[[ch]]
     W[[count]] <- tempW[[ch]] 
     CC[[count]] <- tempCC[[ch]]

     if(length(ch)==1 && tempD[ch]<D[length(D)]) { #compare with old D and store updated values if better
      anyBetter <- TRUE #we will take a new round and check all loci again
     }
     D2 <- tempD[ch]#/length(unlist(Y))
     mhat2 <- tempMx[ch,]
     D <- c(D,D2)
     mhat <- rbind(mhat,mhat2)
     count <- count + 1 #update to next loci-counter
    } #end loci in block: j
   } #end block: i

   show(paste('Mx:',paste(round(mhat2,2),collapse=',',sep=''),'| MD:', round(D2,2))) #show tau:
   if(!anyBetter) break #stop iteration if no loci gave improvement
  } #end repeat{}

  #note: Only one dimension of combinations. So no prioritylist as output
  order <- rep(NA,nL)
  for(i in 1:nL) order[i] <- which(loc==i) #find order of orginal loci
  CC <- matrix(unlist(CC)[order],ncol=nL)  #need correct order of loci
  Mx <-  matrix(mhat[nrow(mhat),],ncol=nC)
  pList <- list()
  pList$MD <- D[length(D)] #select last
  CClist = getCClist(CC,Alist,Mx,ordered=FALSE)
  pList$CC <- CClist$CC
  pList$Mx <- CClist$Mx

  #get simple analysis in addition
  simpleList = mastermix_simple(mixData2,nC,eps=Inf,model,refData2, condOrder, pList=FALSE,zeroMx)$simpleList
  return(list(simpleList=simpleList,pList=pList))
 }  

 #############################################################
 #**********************Peal-off Search *********************#
 #############################################################
 #The algorithm merges loci with decreasing order of number of alleles and minimum peak height and keep the best local subset of combinations (extension of greedy)
  #Note 1: For seperate loci, if more than 10 combinations, the best clustered combinations (distance in MD) are keeped.
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
 mastermix_peeloff = function(mixData2,nC,eps=100,model=1,refData2=NULL, condOrder=NULL,zeroMx=FALSE) {
  #model={1,2,3} - OLS, WLS or GLS
  #mixData: Uses only replicate 1 and contains alleles and heights
  #nC: number assumed contributors at loci
  #refData: list with combinations from reference profiles (in alleles)
  #      -supports a subset of Loci (but both alleles must be in Evidence)
  #condOrder - vector for the specified references. specifies their order.
  #eps: version1:proportion of extracted candidates, version2:max number of combinations
  if(is.null(eps)) stop("eps was not specified!")
  if(eps<=0 | round(eps)!=eps) stop("eps must be a positive integer!")

  cc = 0 #counter of number of loci
  Alist = list() #listed allele info
  Ylist = list() #listed height info
  for(i in 1:length(mixData2$lociname) ) { #for each loci
   cc = cc + 1
   Alist[[cc]] = mixData2$adata[[i]]
   Ylist[[cc]] = mixData2$hdata[[i]]
  }
  nL = cc #number of loci
  nA = unlist(lapply(Ylist,length)) #number of allele
  F= function(x) { if(length(x)==0) { return(0) } else { return(min(x)) } }
  minH = unlist( lapply(Ylist,F) )#minimum allele height
  n = sum(nA) #total number of alleles
  YiSum = unlist(lapply(Ylist,sum))
  Wi_list =  getCovmod(Ylist,model)  #get inverse covariance matricees
 
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
     if(nA[i]<=1) next #should skip if 1 or less allele
     cc = cc + 1
     B[[cc]] <- i
    }
    nbunchs <- nL2 <- cc 
   }
   ##########################
   #fit model for each bunch#
   ##########################
   for(b in 1:nbunchs) { 
    if(seqcount==1) { 
     #Symmetry necessary: Using those with decreasing Mx
     contr_combs <- getContr_combs(Alist,B[[b]],nC,symmetry=TRUE,refs=refData2,condOrder=condOrder)#get combinations for loci
     C[[b]] = list(contr_combs)
    } else if(b>1) {
     next #only modelfit for bunch 1 
    }
    cc <- C[[b]] #is a list itself
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
      Pi <- matrix(getPtab(cc[[i]][j,]),ncol=nC,nrow=nA[B[[b]][i]])
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
     if(j%%1000==0) print(paste('comb',j,'finished for bunch',b))
    } #end each combination
    #all combinations for bunch b calculated:

    #If number of loci in bunch is all
    if(nLinB==nL2) { 
     done=TRUE
     break
    }

    ##############
    ###PEAL OFF###
    ##############
    #peal off obvious errors:
    #hist(MD[keep])
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
    if(sum(keep)==0) { stop('The input of condOrder gave invalid model fit'); }
    if(seqcount==1 && sum(keep)>10) { #if first round: keep lowest cluster (in distance)
     keep2 <- cutree(hclust(dist(MD[keep]),method='single'),2)
     keep2 <- keep2==which.min( c(mean(MD[keep][keep2==1]),mean(MD[keep][keep2==2])) ) 
     keep[keep] <- keep2
    }
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
    #This version (2) ranks loci with number of alleles and minimum height, like the greedy algorithm
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
    print(paste('#Combinations:',paste(newNCombs,collapse=' ')))
    nbunchs <- nB
    B <- newB
    C <- newC
    seqcount <- seqcount + 1
   }
  } #end repeat
  print(paste('Finished! Now postfixing tables.')) 
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

  #get simple analysis in addition
  simpleList = mastermix_simple(mixData2,nC,eps=Inf,model,refData2, condOrder, pList=FALSE,zeroMx)$simpleList
  return(list(simpleList=simpleList,pList=pList))
 } #END PEAL OFF
 

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
    refA = c(refA, refData2$adata[[k]][[i]] )
   }
   #get number of references on locus i:
   nR = sum(condOrder>0)
   if(!is.null(locsel_Ref)) nR = sum( locsel_Ref[i,] &  condOrder>0) #if unselected markers in ref.
   refA = unique(refA) #unique refererence alleles (uncselected are not counted)
   Ai = mixData2$adata[[i]]
   leftA = Ai[!Ai%in%refA] #undescribed alleles for unknown

   if(length(leftA)>2*(nC-nR)) {
    msg <- paste('For locus ',i,', number of unique allele left is ',length(leftA), ', while number of unknowns are ',nC-nR,'. Specify more unknowns or change reference conditioning.',sep='')
    stop(msg)
   }
  }
 }
 if(max(nA)>(2*nC)) {
    msg <- paste('Max alleles in a locus is ',max(nA),'. You must specify a greater number of contributors',sep='')
    stop(msg)
    return(NULL)
 } else {
  if(any(length(grep(tolower(method),'simple')))) f = mastermix_simple
  if(any(length(grep(tolower(method),'greedy'))))  f = mastermix_greedy
  if(any(length(grep(tolower(method),'peeloff')))) f = mastermix_peeloff
  deconvlist <- f(mixData2,nC,eps,model,refData2,condOrder,zeroMx=zeroMx)
 }

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
  ln <- deconvlist$lociname
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
  ln <- deconvlist$lociname
  rownames(combs) <- paste('#',c(t(replicate(nC,1:nJ))),':C',1:nC,sep='')
  ln[nchar(ln)==1] <- paste('0',ln[nchar(ln)==1],sep='')
  colnames(combs) <- c(ln,'Mx','MD') #paste(c(t(replicate(2,ln))),c('_1','_2'),sep='')
  return(combs)
 }

 #add extra result info:
 deconvlist$locinames <- mixData2$lociname #store locinames
 deconvlist$result1 = getCombs(deconvlist)
 deconvlist$result2 = getCombs2(deconvlist)
 deconvlist$data <- list(mixData=mixData2,refData=refData2,locsel_Mix=locsel_Mix,locsel_Ref=locsel_Ref)
 deconvlist$options <- list(condOrder = condOrder,nContr=nC,method=method,model=model,epsilon=eps)
 return(deconvlist)
} 
#***********************END MASTER MIX**************************#

