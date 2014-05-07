###
#TODO

###########
#Changelog#
###########

#05.05.14 - Updated revival of project work.
#30.04.14 - Finished implemented datbase search
#16.04.14 - Bugs: Fix invariant Locus-names,Handle if refData is empty
#21.01.14 - Bugs found before deconvolving
#20.01.14 - Quant-probability function contProb finished implemented (not inserted into GUI).
#25.07.13 - Fix problem with condOrder argument used with keepElite
#14.05.13 - Start create GUI for mastermix-functions (replaces mixinttool)

#' @title mastermixTK
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description mastermixTK is a GUI for the function mastermix for performing linear deconvolution procedure of STR DNA mixtures.
#' @export
#' @usage mastermixTK(envirfile=NULL)
#' @details The function is a graphical layeFr for the function mastermix. See ?mastermix for more information.
#' 
#' The structure of profiles is: Data[['samplename']]$adata[['locusname']] for allele information and Data[['samplename']]$hdata[['locusname']] for allele height information (if any)
#' 
#' @keywords deconvolution, optimization


 #setwd("C:/Users/oebl/Dropbox/Forensic/MixtureProj/myDev/examples")
 #setwd("C:/Users/oebl/Dropbox/Forensic/MixtureProj/myDev/mastermix/R")
 #setwd("D:/Dropbox/Forensic/MixtureProj/myDev/mastermix/R")
 # source("mastermixTK.R");mastermixTK()
 #rm(list=ls()) #must be removed after debugging
 #source("plotEPG.R")
 #source("calcDOdistr.R")
 #source("deconvolve.R")
 #source("getContrCombs.R")

mastermixTK = function(envirfile=NULL) {
 #input: envirfile - a file to saved environment of a project

 #size of main window
 mwH <- 1000
 mwW <- 1500

 #type of gwidgets-kit
 options(guiToolkit="tcltk")

 #Required in deconvolve-function:
 require(gWidgetstcltk) #requires only gWidgets also:

 #version:
 version = 1


 #create own environment within mastermixTK - function (i.e. env-object is enclosed from this function)
 if(is.null(envirfile)) {
  mmTK = new.env() #create new envornment object
  #initializing environment variables
  assign("workdir",NULL,envir=mmTK) #assign working directory to mmTK-environment
  assign("freqfolder",NULL,envir=mmTK) #assign freqfolder to mmTK-environment
  assign("kits",NULL,envir=mmTK) #assign kitList to mmTK-environment
  assign("selPopKitName",NULL,envir=mmTK) #selected kit and population for popFreq
  assign("popFreq",NULL,envir=mmTK) #assign popFreq to mmTK-environment
  assign("mixData",NULL,envir=mmTK) #assign mixdata to mmTK-environment
  assign("refData",NULL,envir=mmTK) #assign refdata to mmTK-environment
  assign("dbData",NULL,envir=mmTK) #assign dbData: matrix referenceDatabase to mmTK-environment (encoded)
  assign("DBsearch",NULL,envir=mmTK) #assign database results
  assign("deconvlist",NULL,envir=mmTK) #assign deconvolved result to mmTK-environment
  assign("LRopt",NULL,envir=mmTK) #assign LR options (weighting)
  assign("LRoptDB",NULL,envir=mmTK) #assign LR options (database)
  assign("freqSize",1000,envir=mmTK) #Assuming N=1000 samples
 } else {
  load(envirfile) #loading environment
 }

 prim = c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113, 127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263, 269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421, 431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593, 599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757, 761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941, 947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093, 1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249, 1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427, 1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549) 

 #minimum frequency
 getminFreq = function(size) {
  return(5/(2*size))
 }

 #Function to get data from environment
 getData = function(type) {
   Data <- NULL
   if(type=="mix") Data <- get("mixData",envir=mmTK) #assign kit to mmTK-environment
   if(type=="ref") Data <- get("refData",envir=mmTK) #assign kit to mmTK-environment 
   if(type=="db") Data <- get("dbData",envir=mmTK) #assign kit to mmTK-environment 
   return(Data)
 }

 #function for inserting mix/ref/db-names into existing gcheckboxgroup
 restoreCheckbox = function(type) {
  subD <- getData(type)
  if(!is.null(subD)) { return( names(subD))
  } else { return(numeric()) }
 }

 #Load kit with allelefreq-info from filefolder:
 loadKitList = function(freqpath) {
  freqfiles = list.files(freqpath)
  kitList <- list()
  for(i in 1:length(freqfiles)) {
   filename = paste(freqpath,"/",freqfiles[i],sep="")
   tab=tableReader(filename)
   Anames = tab[,1]
   tab = tab[,-1]
   freqlist = list()
   for(j in 1:ncol(tab)) { #for each locus
     tmp = tab[,j]
     tmp2 = tmp[!is.na(tmp)]
     names(tmp2) = Anames[!is.na(tmp)]
     freqlist[[j]] = tmp2
   }
   names(freqlist) = toupper(colnames(tab)) #LOCUS-names are assigned as Upper-case
   kit = unlist(strsplit(freqfiles[i],"_"))[1]
   pop = unlist(strsplit(freqfiles[i],"_"))[2]
   pop = unlist(strsplit(pop,"\\."))[1]
   kitind = kit==names(kitList) 
   kitList[[kit]][[pop]] = freqlist
  }
  assign("kits",kitList,envir=mmTK) #assign kit to mmTK-environment
 }

#################
##HELPFUNCTIONS##
#################

 #Function which takes rownames and adds to first column
 addRownameTable = function(tab) {
  tmp <- colnames(tab)
  tab <- cbind(rownames(tab),tab)
  colnames(tab) <- c("X",tmp)
  return(tab)
 }

 #Robust function for reading tables:
 tableReader=function(filename) {
  tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE)
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE)
  return(tab) #need dataframe to keep allele-names correct!!
 }

 #save result table to file:
 saveTable = function(tab,sep="txt") {
  tabfile  = gfile(text="Save table",type="save") #csv is correct format!
  if(!is.na(tabfile)) {
   if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
   if(sep=="txt" | sep=="tab") write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) 
   if(sep=="csv") write.table(tab,file=tabfile,quote=FALSE,sep=",",row.names=FALSE) 
   print(paste("Table saved in ",tabfile,sep=""))
  }
 } #end file

 strsplit2 <- function(x,spl) {
  if(nchar(x)==0) return("")
  txt <- x
  for(j in 1:length(spl)) {
   txt <- unlist(strsplit(txt,split=spl[j]))
  }
  return(txt)
 }

 #Helpfunctions for converting profiles from list to table.
 sample_listToTable = function(Y) {
   sn = names(Y) #Y is a list on form Y[[samplename]]$adata[[locusname]],Y[[samplename]]$hdata[[locusname]]
   aM = 0   #count number of max allele data:
   hM = 0   #count number of max allele heights:
   for(k in 1:length(sn)) { #for each sample
    aM = max(unlist(lapply(Y[[sn[k]]]$adata,length)),aM)
    if(!is.null(Y[[sn[k]]]$hdata)) hM = max(unlist(lapply(Y[[sn[k]]]$hdata,FUN=function(x) { return( sum(x!="")) })),hM)
   }
   #create tables:
   X=numeric()
   for(k in 1:length(sn)) { #for each sample
    newsample=numeric() #for allele
    ln = names(Y[[sn[k]]]$adata)
    for(i in 1:length(ln)) {
     newrow = Y[[sn[k]]]$adata[[ln[i]]]
     newsample = rbind(newsample, c(newrow,rep("",aM-length(newrow))))
    }
    newsample2=numeric() #for heights
    if(hM>0) {
     for(i in 1:length(ln)) {
      newrow = Y[[sn[k]]]$hdata[[ln[i]]]
      newsample2 = rbind(newsample2, c(newrow,rep("",hM-length(newrow))))
     }
    }
    X = rbind(X,cbind(sn[k],ln,newsample,newsample2))
   }
   cn = c("SampleName","Marker", paste("Allele",1:aM,sep=""))
   if(hM>0) cn = c(cn,paste("Height",1:hM,sep=""))
   colnames(X)  = cn
   return(X)
 } #end of functions

 #Helpfunctions for converting profiles from table to list.
 sample_tableToList = function(X) {
  cn = colnames(X) #colnames 
  lind = grep("marker",tolower(cn),fixed=TRUE) #locus col-ind
  sind = grep("sample",tolower(cn),fixed=TRUE) #sample col-ind
  A_ind = grep("allele",tolower(cn),fixed=TRUE) #allele col-ind
  H_ind = grep("height",tolower(cn),fixed=TRUE) #height col-ind
  ln = unique(toupper(X[,lind])) #locus names: Convert to upper case
  sn = unique(X[,sind]) #sample names
  I = length(ln)
  Y = list() #insert non-empty characters:
  for(k in 1:length(sn)) { #for each sample in matrix
   Y[[sn[k]]] = list() #one list for each sample
   if(length(A_ind)>0) Y[[sn[k]]]$adata=list()
   if(length(H_ind)>0) Y[[sn[k]]]$hdata=list() #can also be true for references
   for(i in 1:I) { #for each locus
     xind = X[,sind]==sn[k] & toupper(X[,lind])==ln[i] #get index in X for given sample and locus
     keep <- which(!is.na(X[xind,A_ind]) & X[xind,A_ind]!="")
     if(length(keep)>0) {
      if(length(A_ind)>0) Y[[sn[k]]]$adata[[ln[i]]] = as.character(X[xind,A_ind][keep])
      if(length(H_ind)>0) Y[[sn[k]]]$hdata[[ln[i]]] = as.numeric(as.character(X[xind,H_ind][keep]))
     }
   }
  }
  return(Y)
 }

#Function for database searching for given model in LRsetup:
#should merge database samples with same locus to save calculations!
doDBsearch <- function(LRopt,verbose=TRUE) {
  LRopt$LRfit <- calcLR(LRopt,doLR=FALSE,verbose=FALSE) #get need info of LRcalculation 
  assign("LRoptDB",LRopt,envir=mmTK) #store
#  LRopt <- get("LRoptDB",envir=mmTK)
   maxC <- max(LRopt$uHp+length(LRopt$Hp),LRopt$uHd+length(LRopt$Hd))

  print("Precalculation for database search...")
  LRfit <- LRopt$LRfit #need fittet LR to get ref-data in hypothesis
  locs <- names(LRfit) #get fitted loci (relevant)
  popFreq <-  get("popFreq",envir=mmTK)[locs] #extract relevant loci
  dbData <- getData("db") #get reference database(s)
  dbSel <- names(dbData)

  outD <- nLocs <- macDB <- lrDB <- numeric() #variables to store results of all databases
  for(loc in locs) { #for each locus to consider in mixture
    print(paste("Calculations for locus:",loc))
    Ei <- LRfit[[loc]]$evid #extract sample
    subfreq <- popFreq[[loc]] #extract frequencies
    Ainfo <- names(subfreq) #extract allele-info
    Pinfo <- prim[1:length(Ainfo)]
    subDB <- numeric()
    for(dsel in dbSel) { #for all databases
     subD <- dbData[[dsel]]
     if(loc==locs[1]) {
      outD  <- c(outD,rownames(subD)) #sample names will be character
      M <- nrow(subD) #number of samples in database  
      nLocs <- c(nLocs,rep(0,M))
      macDB <- c(macDB,rep(0,M))
      lrDB <- c(lrDB,rep(0,M))
     }
     dblocs <- toupper(colnames(subD)) #get database locs
     dblocind <- grep(toupper(loc),toupper(dblocs),fixed=TRUE)
     if(all(Ei==0) | length(dblocind)==0 )  next  #check if loci was calculated
     #translate to genotypes
     subDB <- c(subDB,subD[,dblocind]) #merge info of databases
    }
    uniqRef <- unique(subDB) #get unique genotypes
    uniqRef <- uniqRef[!is.na(uniqRef)] #remove NA

    for(j in 1:length(uniqRef)) { #for each unique genotypes      
     rowind <- which(subDB==uniqRef[j]) #samples with this genotype
     DBreff <- Ainfo[uniqRef[j]%%Pinfo==0] #get allele-info
     if(length(DBreff)==1) DBreff <- rep(DBreff,2) #homozygote genotype
     hpref = c(DBreff,LRfit[[loc]]$refHp) #reff under hp
     hdref = LRfit[[loc]]$refHd #reff under hd
     #quick calculation if Q-assignation
     if(LRopt$Qcalc) {
      subfreq <- subfreq[names(subfreq)%in%unique(c(Ei,hpref,hdref))]
      subfreq <- c(subfreq ,1-sum(subfreq))
     }
     hp = likEvid(Repliste=Ei,T=hpref,V=NULL,x=LRopt$uHp,theta=LRopt$theta,prDHet=rep(LRopt$DOprob,maxC),prDHom=rep(LRopt$DOprob^2,maxC),prC=rep(LRopt$DIprob,maxC),freq=subfreq)
     hd = likEvid(Repliste=Ei,T=hdref,V=DBreff,x=LRopt$uHd,theta=LRopt$theta,prDHet=rep(LRopt$DOprob,maxC),prDHom=rep(LRopt$DOprob^2,maxC),prC=rep(LRopt$DIprob,maxC),freq=subfreq)
     lrDB[rowind] <- lrDB[rowind] + log10(hp/hd)
     macDB[rowind] <- macDB[rowind] + sum(DBreff%in%Ei) #count number of fitting evidence in compound evidence
     nLocs[rowind] <- nLocs[rowind] + 1
    } #end for each unique ref
   } #end for each locus
  ord <- order(lrDB,decreasing=TRUE)
  result <- cbind(outD[ord],lrDB[ord],macDB[ord],nLocs[ord]) #extend
  colnames(result) <- c("Reference","log10LR","MAC","nLocs")
  assign("DBsearch",result,envir=mmTK)
#  dbwin <- gwindow("Result of database search", visible=FALSE, width=mwW,height=mwH)
#  gtable(outD,container=dbwin,multiple=TRUE) #create table
#  visible(dbwin) <- TRUE
}

#this is only used for tippet!
getLRRMlist = function(LRopt,tippetSel) {
 print(paste("Tippet Man selected:",tippetSel))
# LRopt <-  get("LRopt",envir=mmTK) #extract relevant loci
 LRfit <- LRopt$LRfit #need fittet LR to get ref-data in hypothesis
 locs <- names(LRfit) #get fitted loci
 popFreq <-  get("popFreq",envir=mmTK)[locs] #extract relevant loci
 #prepare random profile generator:
 LRRM <- list()
 Glist <- list()
 H_p = list()
 H_d = list()
 for(loc in locs) {
  Ei <- LRfit[[loc]]$evid #extract sample
  if(all(Ei==0))  next  #check if loci was calculated
  print(paste("RMLR calculation for locus:",loc,"..."))
  maxC <- max(LRopt$uHp+length(LRopt$Hp),LRopt$uHd+length(LRopt$Hd))
  subfreq <- popFreq[[loc]]

  G = t(as.matrix(expand.grid(rep(list(as.numeric(names(subfreq)),as.numeric(names(subfreq)) )))))
  keep = G[2,]>=G[1,] #unique genotypes
  G <- G[,keep]  #store genotypes
  G <- matrix(as.character(G),nrow=2) #make string names
  tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(subfreq),as.numeric(subfreq) )))))
  Gprob = exp(colSums(log(tmpP[,keep]))) #get allele probs
  ishet = G[1,]!=G[2,]
  Gprob[ishet] = 2*Gprob[ishet] #multiply with two to get heterozygote prob

  Glist[[loc]] <- list(G=G,Gprob=Gprob)
  H_p[[loc]] = rep(NA,ncol(G))
  H_d[[loc]] = rep(NA,ncol(G)) #must be considered!
  OindHp <- which(!colnames(LRfit[[loc]]$refHp)%in%tippetSel) #get references not for tippet 
  OrefsHp <- c(LRfit[[loc]]$refHp[,OindHp])  #get other alleles than tippet
  OrefsHd <- LRfit[[loc]]$refHd 
  for(j in 1:ncol(G)) { #for each genotype
   reff = c(G[,j],OrefsHp)
   if(LRopt$Qcalc) { #quick calculation if Q-assignation
    subfreq <- subfreq[names(subfreq)%in%unique(c(Ei,reff))]
    subfreq <- c(subfreq ,1-sum(subfreq))
   }
   H_p[[loc]][j] = likEvid(Repliste=LRfit[[loc]]$evid,T=reff,V=NULL,x=LRopt$uHp,theta=LRopt$theta,prDHet=rep(LRopt$DOprob,maxC),prDHom=rep(LRopt$DOprob^2,maxC),prC=rep(LRopt$DIprob,maxC), freq=subfreq)
   H_d[[loc]][j] = likEvid(Repliste=LRfit[[loc]]$evid,T=OrefsHd,V=G[,j],x=LRopt$uHd,theta=LRopt$theta,prDHet=rep(LRopt$DOprob,maxC),prDHom=rep(LRopt$DOprob^2,maxC),prC=rep(LRopt$DIprob,maxC), freq=subfreq)
  } #end for:j genotypes
  LRRM[[loc]] <- H_p[[loc]]/H_d[[loc]] #constant value
 } #end for each loci
 print("..RMLR precalculation done!")
 return(list(LRRM=LRRM,Glist=Glist,H_p=H_p,H_d=H_d))
} 

#Helpfunctions for calculating p-value
calcPvalue = function(LRRMlist,lrobs) {
 #LRRMlist - object returned from getLRRMlist
 #lrobs - observed LR to find p-value of
 lrlist = LRRMlist$LRRM #get list with LR-calculations
 K <- length(lrlist) #number of loci
 M <- max(sapply(lrlist,length)) #max number of RM
 X <- matrix(0,K,M)
 P <- matrix(0,K,M)
 cc = 0 #counter for matrices
 for(i in 1:K) {
  if(is.null(lrlist[[i]])) next
  cc = cc + 1
  idx <- order(lrlist[[i]],decreasing=TRUE)
  X[cc,1:length(lrlist[[i]])] <- lrlist[[i]][idx]
  P[cc,1:length(lrlist[[i]])] <- LRRMlist$Glist[[i]]$Gprob[idx]
 } 
 X = X[1:cc,]
 P = P[1:cc,]
 print("Calculating p-value...")
 outC <- pvalueC(lrobs,X,P) 
 return(outC) #return P-value
}


###################################################################
###########################GUI#####################################
###################################################################

 #print LR model in guitab-layout:
   printModel = function(LRopt,guitab) { #function for plotting function
    nM <- length(LRopt$Evidsel)
    nonHdcon <- LRopt$Hp[!LRopt$Hp%in%LRopt$Hd] #get possible tippets
    leftHp <- LRopt$Hp[LRopt$Hp%in%LRopt$Hd] #get possible tippets
    nNotHd <- length(nonHdcon) #number of possible tippets
    nHp <- length(LRopt$Hp)
    nleftHp <- length(leftHp)#nHp-nNotHd #non-tippets
    nHd <- length(LRopt$Hd)
    guitab[1,1] <- glabel(text="Evidence(s):",container=guitab)
    for(nm in 1:nM) {
     guitab[1+nm,1] <- glabel(text=names(LRopt$Evidsel)[nm],container=guitab)
    }
    guitab[2+nM,1] <- glabel(text="",container=guitab) #space
    guitab[3+nM,1] <- glabel(text="Contributor(s) under Hp:",container=guitab)
    if(nNotHd>0) {
      guitab[4+nM,1] <- gradio(items=nonHdcon,container=guitab,checked=FALSE)
    }
    if(nleftHp>0) {
      guitab[5+nM,1] <- gradio(items=leftHp,container=guitab,checked=FALSE)
      enabled(guitab[5+nM,1]) <- FALSE
    }
    guitab[6+nM,1] <- glabel(text="",container=guitab) #space
    guitab[7+nM,1] <- glabel(text="Contributor(s) under Hd:",container=guitab)
    if(nHd>0) {
     guitab[8+nM,1] <- gradio(items=LRopt$Hd,container=guitab,checked=FALSE)
     enabled(guitab[8+nM,1]) <- FALSE
    }
    guitab[9+nM,1] <- glabel(text="",container=guitab) #space
    guitab[10+nM,1] <- glabel(text="",container=guitab) #space
   
    #Parameter selection:
    guitab[11+nM,1] <- glabel(text="Parameters:",container=guitab)
    guitab[12+nM,1] <- glabel(text=paste("#unknowns (Hp):",LRopt$uHp),container=guitab)
    guitab[13+nM,1] <- glabel(text=paste("#unknowns (Hd):",LRopt$uHd),container=guitab)
    guitab[14+nM,1] <- glabel(text=paste("Drop-out prob: ",LRopt$DOprob),container=guitab)
    guitab[15+nM,1] <- glabel(text=paste("Drop-in prob: ",LRopt$DIprob),container=guitab)
    guitab[16+nM,1] <- glabel(text=paste("Theta:",LRopt$theta),container=guitab)
    guitab[17+nM,1] <- glabel(text=paste("Q-assignation:",LRopt$Qcalc),container=guitab)
   }

 #Menu bar file-lists:
 f_setwd = function(h,...) {
  dirfile = gfile(text="Select folder",type="selectdir")
  if(!is.na(dirfile)) {
   setwd(dirfile)
   assign("workdir",dirfile,envir=mmTK) #assign working directory
  }
 }
 f_openproj = function(h,...) {
  projfile = gfile(text="Open project",type="open", filter=list("Project"=list(patterns=list("*.Rdata"))))
  if(!is.na(projfile)) {
   dispose(mainwin)
   mastermixTK(projfile) #send environment into program
  }
	 }
 f_saveproj = function(h,...) {
  projfile = gfile(text="Save project",type="save")
  if(!is.na(projfile)) {
   if(length(unlist(strsplit(projfile,"\\.")))==1) projfile = paste0(projfile,".Rdata")
   print(sapply(mmTK,object.size)/1e6)
   save(mmTK,file=projfile) #save environment
   print(paste("Project saved in ",projfile,sep=""))
  }
 }
 f_quitproj = function(h,...) {
  ubool <- gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(svalue(ubool)) {
   savefile <- gfile(text="Save file",type="save")
   if(length(unlist(strsplit(savefile ,"\\.")))==1) savefile = paste0(savefile ,".Rdata")
   save(mmTK,file=savefile) 
   print(paste("Project saved in ",savefile,sep=""))
  } else { 
   print("Program terminated without saving")
  }
  dispose(mainwin) #remove window
 }

###########GUI WINDOW STARTS#########

 #Menu bar - tags
 mblst = list( #project saving and so on
  File=list(  
   'Set directory'=list(handler=f_setwd),
   'Open project'=list(handler=f_openproj),
   'Save project'=list(handler=f_saveproj),
   'Quit'=list(handler=f_quitproj,icon="close")
  )
 )
 
##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 #change working directory to the one stored in mmTK-environment
 wd=get("workdir",envir=mmTK) #assign working directory to mmTK-environment
 if(!is.null(wd)) setwd(wd)
 
 #Main window:
 spc <- 10
 mainwin <- gwindow(paste("Mixture analysis Tool v",version,sep=""), visible=FALSE, width=mwW,height=mwH)
 gmenu(mblst,container=mainwin)
 nb = gnotebook(container=mainwin)
 tab1 = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Create data")
 tab2 = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Import data")
 tab3 = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Deconvolution Setup")
 tab4 = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Deconvolution Results")
 tab5 = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="LR Setup")
 tab6 = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="LR Results")
 tab7 = ggroup(horizontal=FALSE,expand=TRUE,spacing=spc,container=nb,label="Database Search")
 svalue(nb) <- 2 #initial start at second tab

####################################################
###############Tab 1: Create Data:##################
####################################################
 #layout: 
 tab1a = glayout(spacing=3,container=gframe("Generator",container=tab1))
 tab1b = glayout(spacing=3,container=gframe("Edit",container=tab1)) #create dataframe
 tab1c = glayout(spacing=20,container=gframe("Import/Export profile",container=tab1)) #store dataframe

 refreshTab1 = function() {
  #Function for achiving entry-text: Reads into list-format:
  getEntryDataGUI <- function(asList=TRUE) {
   sname <- svalue(tab1c[1,2])
   Data <- list()
   Data[[sname]] = list(adata=list(),hdata = list())
   ln = numeric()
   for(i in 1:(length(get("popFreq",envir=mmTK))+1)) { #include AMEL in last entry
     ln <- c(ln,svalue(tab1b[i+1,1]))
     Data[[sname]]$adata[[i]] <- strsplit2(svalue(tab1b[i+1,2]),c(" ",",",";",":"))
     Data[[sname]]$hdata[[i]] <- strsplit2(svalue(tab1b[i+1,3]),c(" ",",",";",":"))
   }
   names(Data[[sname]]$adata) = ln
   names(Data[[sname]]$hdata) = ln
   return(Data)
  }
  #Function for achiving entry-text: Write from list-format:
  setEntryDataGUI <- function(Data,sname) {
   svalue(tab1c[1,2]) <- sname #insert samplename to GUI
   for(i in 1:length(Data[[sname]]$adata)) {
    svalue(tab1b[i+1,2]) <- paste(Data[[sname]]$adata[[i]],collapse=",",sep="")
    if(is.null(Data[[sname]]$hdata)) {
     svalue(tab1b[i+1,3]) <- ""
    } else {   
     svalue(tab1b[i+1,3]) <- paste(Data[[sname]]$hdata[[i]],collapse=",",sep="")
    }
   }
  }
  #Randomize profile: Simulate random and set entryData
  #could set type: Reference or Mixtur from h$type
  ranProfile = function(h,...) {
   freqList <- get("popFreq",envir=mmTK)  #get selected population frequency
   sname <- "random" #samplename
   Data <- list(list(adata=list()))
   names(Data) <- sname
   if(length(freqList)>0) {
    for(i in 1:length(freqList)) {  
     G = t(as.matrix(expand.grid(rep(list(as.numeric(names(freqList[[i]])),as.numeric(names(freqList[[i]])) )))))
     keep = G[2,]>=G[1,] #unique genotypes
     G <- G[,keep]  #store genotypes
     tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(freqList[[i]]),as.numeric(freqList[[i]]) )))))
     Gprob = exp(colSums(log(tmpP[,keep]))) #get allele probs
     ishet = G[1,]!=G[2,]
     Gprob[ishet] = 2*Gprob[ishet] #multiply with two to get heterozygote prob

     if(h$action=="ref") {
      Data[[sname]]$adata[[names(freqList)[i]]] <- G[,sample(1:length(Gprob),size=1,prob=Gprob,replace=TRUE)]
     }
     if(h$action=="mix") {
      sampledA = c(G[,sample(1:length(Gprob),size=2,prob=Gprob,replace=TRUE)])
      sampledH = round(abs(rnorm(length(sampledA),1000,100)))
      agg <- aggregate(sampledH,by=list(sampledA),sum)
      Data[[sname]]$adata[[names(freqList)[i]]] <- agg[,1]
      Data[[sname]]$hdata[[names(freqList)[i]]] <- agg[,2]
     }     
    } 
   }
   if(h$action=="ref")  amsample = sample(1:2,size=1,replace=TRUE)
   if(h$action=="mix")  amsample = sample(1:2,size=2,replace=TRUE)
   ams <- as.numeric()
   for(k in 1:length(amsample)) {
     if(amsample[k]==1) ams= c(ams,"X","Y")
     if(amsample[k]==2) ams= c(ams,"X","X") 
   }
   if(h$action=="ref")  Data[[sname]]$adata[["AMEL"]] <- ams
   if(h$action=="mix") {
      sampledH = round(abs(rnorm(length(ams),1000,100)))
      agg <- aggregate(sampledH,by=list(ams),sum)
      Data[[sname]]$adata[["AMEL"]] <- agg[,1]
      Data[[sname]]$hdata[["AMEL"]] <- agg[,2]
   }
   print(Data)
   setEntryDataGUI(Data,sname)  
  } #end

  #load/save profiles
  f_openprof = function(h,...) {
   # proffile = gfile(text=paste("Open profile",sep=""),type="open")
    proffile = gfile(text="Open profile",type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab"))))
    if(!is.na(proffile)) {
     Data = tableReader(proffile) #load profile
     print(Data)
     Data = sample_tableToList(Data) #convert from table to list
     setEntryDataGUI(Data,names(Data)[1]) #set imported profile to GUI
    }
  }
  f_saveprof = function(h,...) {
    Data = getEntryDataGUI(asList=FALSE) #get profile data as table
    Data = sample_listToTable(Data) #convert from table to list
    saveTable(Data,"csv")
  }

  #Layout2
  tab1a[1,1] = gbutton(text="Generate ref-profile",container=tab1a,handler=ranProfile,action="ref")
  tab1a[1,2] = gbutton(text="Generate 2-mixed-profile",container=tab1a,handler=ranProfile,action="mix")

  freqList <- get("popFreq",envir=mmTK)  #get selected population frequency
  locnames = c(names(freqList),"AMEL") #add also AMEL
  nL <- length(locnames) #number of loci in kit
  tab1b[1,1] <- glabel("Marker",container=tab1b)
  tab1b[1,2] <- glabel("Allele",container=tab1b)
  tab1b[1,3] <- glabel("Height",container=tab1b)
  for(i in 1:nL) {
   tab1b[i+1,1] <-  glabel(locnames[i],container=tab1b)
   tab1b[i+1,2] <-  gedit("",container=tab1b)
   tab1b[i+1,3] <-  gedit("",container=tab1b)
  }

  tab1c[1,1] = glabel(text="Profile name",container=tab1c)
  tab1c[1,2] = gedit("",container=tab1c)
  tab1c[2,1] = gbutton(text="Open profile",container=tab1c,handler=f_openprof)
  tab1c[2,2] = gbutton(text="Save profile",container=tab1c,handler=f_saveprof)
 } #end refreshTab1

####################################################
###############Tab 2: Import Data:##################
####################################################

 #a) button for loading kits from directory:
 f_loadkd = function(h,...) { loadKitList(freqpath=get("freqfolder",envir=mmTK)); }

 #b) load/save profiles/database: Supports any filetype
 f_importprof = function(h,...) {
  type=h$action #get type of profile
#  proffile = gfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab"))))
  proffile = gfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*."))))
  if(!is.na(proffile)) {
   Data = tableReader(proffile) #load profile
   print("Raw fil import")
   print(Data[1:min(nrow(Data),100),]) #print raw-input data

###################
##DATABASE IMPORT##
###################
   if(type=="db") { #importing database file
    popFreq <- get("popFreq",envir=mmTK) 
    if(is.null(popFreq)) {
     gmessage("Population frequencies needs to be imported for database search",title="Error")
    } else {
     minFreq <- getminFreq(get("freqSize",envir=mmTK)) #get assigned minFrequency
     #saving MEMORY by convert database here!
     #Data <- dbToGenoTable(Data) #convert database table to genotype matrix

     #Assumption: No reference has same samplename. All samplenames are rowed sequentially
     cn = colnames(Data) #colnames 
     lind = grep("marker",tolower(cn),fixed=TRUE) #locus col-ind
     sind = grep("sample",tolower(cn),fixed=TRUE) #sample col-ind
     aind = grep("allele",tolower(cn),fixed=TRUE) #allele col-ind
     aind <- aind[1:2] #assume only 2 possible alles in reference profiles
     locsDB <- unique(Data[,lind]) #unique markers in DB
     locsPop <- toupper(names(popFreq)) #markers in population
     sn <- unique(Data[,sind]) #unique samples in DB
     newData <- numeric() #encoded reference table
     dbDatalocs <- numeric() #locus in newData
     for(i in 1:length(locsDB)) { #for each marker in DB:
      locind <- grep(toupper(locsDB[i]),toupper(names(popFreq)),fixed=TRUE) #get position in popFreq
      if(length(locind)==0) next #if locus not existing in popFreq
      newCol <- rep(NA,length(sn)) #new column in newData
      subA <- Data[Data[,lind]==locsDB[i],aind] #extract allele data with matching marker
      subS <- Data[Data[,lind]==locsDB[i],sind] #extract sample names with matching marker
      isHom <- which(subA[,2]=="" | is.na(subA[,2])) #if homozygote assignation:
      if(length(isHom)>0) subA[isHom,2] <- subA[isHom,1] #copy first allele
      okSind <- which(sn%in%subS) #samples to consider for particular marker
      if(length(okSind)==0) next #if no samples exists
      newCol[okSind] <- 1 #init existing as 1. NA for missing allele-info
      doneEncode <- matrix(FALSE,ncol=2,nrow=length(okSind)) #matrix checking which we are finished with
      Afreqs <- names(popFreq[[locind]]) #get allele-names. Update for each column
      for(k in 1:length(aind)) { #for each allele-column
       for(j in 1:length(Afreqs)) { #for each unique allele in popFreq:
        okAind <- which(subA[,k]==Afreqs[j]) #find matching alleles in subA[okSind]
        if(length(okAind)==0) next
        doneEncode[okAind,k] = TRUE #check that is finished
        newCol[okSind][okAind] <- newCol[okSind][okAind]*prim[j] #multiply with primenumber
       } #end for each allele j

       #THREAT NEW ALLELES,MISSTYPOS ETC:
       if(any(!doneEncode[,k])) { #if any not encoded
        newA <- unique(subA[!doneEncode[,k],k]) #get new alleles
        newA <- newA[!is.na(newA)] #important to remove NA's
        if(length(newA)==0) next
        tmp <- popFreq[[locind]]
        popFreq[[locind]] <- c(tmp, rep(minFreq,length(newA)))
        names(popFreq[[locind]]) <- c(names(tmp),newA) #add unique
        print(paste0("In locus ",locsDB[i],": Allele(s) ",newA," was inserted with min. frequency ",prettyNum(minFreq)))
        for(j in 1:length(newA)) { #for each unique allele in popFreq:
         okAind <- which(subA[,k]==newA[j]) #find matching alleles in subA[okSind]
         if(length(okAind)==0) next
         newCol[okSind][okAind] <- newCol[okSind][okAind] * prim[ which(names(popFreq[[locind]])==newA[j]) ] #multiply with EXTENDED primenumber
        } #end for each allele j
       } #end if not encoded 
      } #end for each column k
      dbDatalocs <- c(dbDatalocs,toupper(names(popFreq)[locind])) #all locus
      newData <- cbind(newData,newCol) #add column
     } #end for each locus i
   
     #RESCALE popFreq?
     for(i in 1:length(popFreq)) {
      popFreq[[i]] <- popFreq[[i]]/sum(popFreq[[i]])
     }
     print("popFreq was normalized")
     colnames(newData) <- dbDatalocs
     rownames(newData) <- sn #insert sample names
     print(paste0("Database successfully imported with ",nrow(newData)," samples."))

     tmp <- unlist(strsplit(proffile,"/",fixed=TRUE)) #just label the file
     fname <- tmp[length(tmp)] #get filename
     fname <- unlist(strsplit(fname,"\\."))[1]

     dbData <- get("dbData",envir=mmTK) #get already existing databases
     if(is.null(dbData)) dbData <- list() #create list-object
     dbData[[fname]] <- newData #add database to list

     assign("dbData",dbData,envir=mmTK) #store matrix in environment for later use
     assign("popFreq",popFreq,envir=mmTK) #assign updated popFreq
     tab2b[2,3][] <- names(dbData)
#     tab2b[2,3][] <- c(tab2b[2,3][],fname)  #if applying existing
#     tab2b[2,3][] <- tmp[length(tmp)] #use latest only
#     enabled(tab2b[2,3]) <- FALSE #all is considered in the database search
    } #end if popFreq exist
   } else { 
    Data = sample_tableToList(Data) #convert from table to list 
    #get already stored data:
    if(type=="mix") Data2 <- getData("mix") #get data from mmTK-environment
    if(type=="ref") Data2 <- getData("ref") #get data from mmTK-environment

    if(is.null(Data2)) { #if no previous already there
     Data2 <- Data
    } else {
     for(k in 1:length(Data)) { #for each profile
      Data2[[names(Data)[k]]] <- Data[[k]] #insert dataframe
     }
    }
    if(type=="mix")  assign("mixData",Data2,envir=mmTK) #assign data to mmTK-environment
    if(type=="ref")  assign("refData",Data2,envir=mmTK) #assign data to mmTK-environment
    if(type=="mix")  tab2b[2,1][] <- names(Data2)
    if(type=="ref")  tab2b[2,2][] <- names(Data2)
   }
  }
 }

 #c) 
 #prints profiles and EPG or a viewData-table
 f_viewdata = function(h,...) {
   selD <- getData(h$action) #get selected data
   if(h$action=="mix") { #call on epg-function for each mixtures
    kit <- svalue(tab2a[2,1]) #get kit name. Must be same name as in generateEPG
    mixSel <- svalue(tab2b[2,1])  #get selected mixtures
    for(msel in mixSel) {
     subD <- selD[[msel]]
     print("------------------------------------")
     print(paste("Samplename: ",msel,sep=""))
     print(subD)
#     for(loc in names(subD$adata)) {
      #if(length(grep("AMEL",loc,fixed=TRUE))>0) next
#      subD$adata[[loc]] <- as.numeric(subD$adata[[loc]])
#     }
     if(msel!=mixSel[1]) dev.new() #create new plot for next EPG
     plotEPG(subD,kit,msel) #plot epg's
    } 
   }
   if(h$action=="ref") { #print tables only
    refSel <- svalue(tab2b[2,2])  #get selected references
    for(rsel in refSel) {
     print("------------------------------------")
     print(paste("Samplename: ",rsel,sep=""))
     print(selD[[rsel]])
    }
   }
   if(h$action=="db") {  #view imported databases
    dbSel <- svalue(tab2b[2,3])  #get selected database
    popFreq <- get("popFreq",envir=mmTK)
    for(dsel in dbSel) { #for each selected databases
     subD <- selD[[dsel]]
     dblocs <- toupper(colnames(subD)) #get database locs
     outD <- rownames(subD) #will be character
     for(i in 1:length(dblocs)) { #for each locus     
      Ainfo <- names(unlist(popFreq[[dblocs[i]]])) #extract allele-info
      #translate to genotypes
      Pinfo <- prim[1:length(Ainfo)]
      G = t(as.matrix(expand.grid(rep(list(Ainfo,Ainfo )))))
      GP = t(as.matrix(expand.grid(rep(list(Pinfo,Pinfo )))))
      keep = GP[2,]>=GP[1,] #unique genotypes
      G <- G[,keep]  #store genotypes
      G <- paste0(G[1,],paste0("/",G[2,]))
      GP <- GP[,keep]  #store genotypes
      GP <- GP[1,]*GP[2,]
      newRow <- rep(NA,nrow(subD)) 
      for(j in 1:length(GP)) { #for each genotype
       rowind <- which(subD[,i]==GP[j]) #samples with this genotype
       if(length(rowind)==0) next
       newRow[rowind] <- G[j]
      } #end for each allele
      outD <- cbind(outD,newRow) #extend
     }
     colnames(outD) <- c("Reference",dblocs)
     if(nrow(outD)>10000) {
      #gmessage(message="Database contained more than 10000 samples.\n Showing first 10000 samples only!")
      print("Database contained more than 10000 samples. Showing first 10000 samples only!")
      outD <- outD[1:10000,]
     }
     dbwin <- gwindow(paste0("References in imported database ",dsel), visible=FALSE,height=mwH)
     gtable(outD,container=dbwin,multiple=TRUE) #create table
     visible(dbwin) <- TRUE
    } #end for each databases
   } #end if db -case
 }  #end viewdata

 ###############
 #start layout:#
 ###############
 tab2a = glayout(spacing=5,container=gframe("Population frequencies",container=tab2)) #kit and population selecter
 tab2b = glayout(spacing=5,container=gframe("Evidence, Reference, Database",container=tab2)) #evidence,ref dataframe
 tab2d = glayout(spacing=20,container=gframe("Interpretations",container=tab2)) #Tasks button

 #Choose box and import button
 tab2a[1,1] = gbutton(text="1) Select directory",container=tab2a,handler=
  function(h,...) {
   dirfile = gfile(text="Select folder",type="selectdir")
   if(!is.na(dirfile)) assign("freqfolder",dirfile,envir=mmTK) #assign freqfolder
  }
 )
 tab2a[1,2] = gbutton(text="2) Import from directory",container=tab2a,handler=
  function(h,...) {
   loadKitList(freqpath=get("freqfolder",envir=mmTK))
   kitList <- get("kits",envir=mmTK)
   tab2a[2,1][] <- names(kitList)
  }
 ) 
 tab2a[2,1] <- gcombobox(items="", width=100, selected =0 , editable = FALSE, container = tab2a, handler=
    function(h,...) {
     if(!is.null(get("dbData",envir=mmTK))) {
      print("You can not change selected kit after loading a database!") 
     } else {
      kitList <- get("kits",envir=mmTK)
      if(!is.null(kitList)) tab2a[2,2][] <- names(kitList[[svalue(tab2a[2,1])]])
     }
    })
 #population-selection
 tab2a[2,2] <- gcombobox(items="", width=100, selected = 0 , editable = FALSE , container = tab2a, handler=
    function(h,...) {
     if(!is.null(get("dbData",envir=mmTK))) {
      print("You can not change selected population after loading a database!") 
     } else {
      kitList <- get("kits",envir=mmTK)
      if(!is.null(kitList)) {
       popList <- kitList[[svalue(tab2a[2,1])]][[svalue(tab2a[2,2])]] #get selected frequencies
       assign("popFreq",popList,envir=mmTK) #assign popFreq get("popFreq",envir=mmTK)
       popkitname <- c(svalue(tab2a[2,1]),svalue(tab2a[2,2])) #store selected kit and popnames
       assign("selPopKitName",popkitname,envir=mmTK) 
      }
     }
    })
#previous stored kit/pop-names:
 popkitname <- get("selPopKitName",envir=mmTK) 
 if(!is.null(popkitname)) {
  tab2a[2,1][] <- popkitname[1]
  tab2a[2,2][] <- popkitname[2]
 }

 #Select number of samples in popFreq-data
 tab2a[1,3] <-  gbutton(text="Update #Samples",container=tab2a,handler=
  function(h,...) {
   freqsize <- as.numeric(svalue(tab2a[2,3])) #updated size
   assign("freqsize ",freqsize ,envir=mmTK) #assign minFrequency. Assuming N=1000 samples
   minFreq <- getminFreq(freqsize) #get new min-frequencies
   print(paste("New alleles imputed with minimum frequency",minFreq))
  })

 defsize = 1000
 freqSize <- get("freqSize",envir=mmTK) #get assigned minFrequency
 if(!is.null(freqSize)) defsize <- freqSize #size previously stored
 tab2a[2,3] <-  gedit(defsize,container=tab2a,width=8) #default is 1000

 #view popFreq-data in a new window
 tab2a[1,4] <-  gbutton(text="View frequencies",container=tab2a,handler=
  function(h,...) {
   popFreq <- get("popFreq",envir=mmTK) #get frequencies
   if(is.null(popFreq)) {
    tkmessageBox(message="Please import and select population frequencies!")
   } else {
    locs <- names(popFreq)
    unAchr <- unique(unlist(sapply( popFreq,names) )) #get character alleles
    ord <- order(as.numeric(unAchr)) 
    unAchr <- unAchr[ord]  #make increasing order
    outD <- unAchr
    for(i in 1:length(locs)) {
     newRow <- rep(NA,length(unAchr))
     for(j in 1:length(popFreq[[i]])) {
      rowind <- which(unAchr==names( popFreq[[i]][j] ))
      newRow[rowind] <- popFreq[[i]][j]
     }
     outD <- cbind(outD,newRow)
    }
    colnames(outD) = c("Allele",locs) 
    dbwin <- gwindow("Population frequencies", visible=FALSE,height=mwH)
    gtable(outD ,container=dbwin,multiple=TRUE) #create table
    visible(dbwin) <- TRUE
   }
  })

 

 #Choose box and import button
 tab2b[1,1] = gbutton(text="Import evidence",container=tab2b,handler=f_importprof,action="mix")
 tab2b[2,1] = gcheckboxgroup(items="", container = tab2b)
 tab2b[2,1][] = restoreCheckbox("mix")

 #Choose box and import button
 tab2b[1,2] = gbutton(text="Import reference",container=tab2b,handler=f_importprof,action="ref")
 tab2b[2,2] = gcheckboxgroup(items="", container = tab2b)
 tab2b[2,2][] = restoreCheckbox("ref")

 #Choose box and import button
 tab2b[1,3] = gbutton(text="Import database",container=tab2b,handler=f_importprof,action="db")
 tab2b[2,3] = gcheckboxgroup(items="", container = tab2b)
 tab2b[2,3][] = restoreCheckbox("db")

 #view data:
 tab2b[3,1] = gbutton(text="View evidence",container=tab2b,handler=f_viewdata,action="mix")
 tab2b[3,2] = gbutton(text="View references",container=tab2b,handler=f_viewdata,action="ref")
 tab2b[3,3] = gbutton(text="View database",container=tab2b,handler=f_viewdata,action="db")

 #Button-choices further:
 tab2d[1,1] = gbutton(text="Create profiles",container=tab2d,handler=
  function(h,...) {
   refreshTab1() #refresh table for creating profiles
   svalue(nb) <- 1 #change tab of notebook
  }
 )
 tab2d[1,2] = gbutton(text="Deconvolution",container=tab2d,handler=
  function(h,...) { 
   mixSel <- refSel <- numeric()
   if(length(tab2b[2,1][])>0) mixSel <- svalue(tab2b[2,1])  #get selected mixtures
   if(length(tab2b[2,2][])>0) refSel <- svalue(tab2b[2,2])  #get selected references
   if(length(mixSel)==0) {
    tkmessageBox(message="Please import and select mixture-profile!")
   } else {
    refreshTab3(mixSel,refSel) #refresh table with selected data
    svalue(nb) <- 3 #change tab of notebook
   }
  }
 ) 
 tab2d[1,3] = gbutton(text="LR calculation",container=tab2d,handler=
  function(h,...) { 
   popFreq <- get("popFreq",envir=mmTK)
   mixSel <- refSel <- numeric()
   if(length(tab2b[2,1][])>0) mixSel <- svalue(tab2b[2,1])  #get selected mixtures
   if(length(tab2b[2,2][])>0) refSel <- svalue(tab2b[2,2])  #get selected references
   if(length(mixSel)==0) {
    gmessage(message="Please import and select mixture-profile!")
   } else if(length(refSel)==0) {
    gmessage(message="Please import and select reference-profile!")
   } else if(is.null(popFreq)) {
    gmessage("No frequencies was specified!\n Please import table.")
   } else {
    refreshTab5(mixSel,refSel,"weight") #refresh table with selected data
    svalue(nb) <- 5 #change tab of notebook
   }
  }
 ) 
 tab2d[1,4] = gbutton(text="Database search",container=tab2d,handler=
  function(h,...) { 
   mixSel <- refSel <- numeric()
   if(length(tab2b[2,1][])>0) mixSel <- svalue(tab2b[2,1])  #get selected mixtures
   if(length(tab2b[2,2][])>0) refSel <- svalue(tab2b[2,2])  #get selected references
   if(length(mixSel)==0) {
    tkmessageBox(message="Please import and select mix-profile!")
   } else if(length(tab2b[2,3][])==0) {
    tkmessageBox(message="Please import database!")    
   } else {
    refreshTab5(mixSel,refSel,"dbsearch") #refresh table with selected data
    svalue(nb) <- 5 #change tab of notebook
   }
  }
 ) 
# enabled(tab2d[1,4]) <- FALSE



############################################################
###############Tab 3: Deconvolution setup:##################
############################################################

 #layout: 
  tab3tmp <- glayout(spacing=30,container=tab3)

 refreshTab3 = function(mixSel,refSel) {
   tab3c = glayout(spacing=1,container=(tab3tmp[1,1] <-gframe("Configurations",container=tab3tmp))) 
   tab3b = glayout(spacing=1,container=(tab3tmp[1,2] <-gframe("Data selection",container=tab3tmp))) 

   #mixSel,refSel is name of profiles in mixData,refData
   #note: User includes those who are interest to select for conditioning directly!
   #check if already exists:
   nM = length(mixSel) #number of mix-profiles
   nR = length(refSel) #number of ref-profiles
   mixD = getData("mix")
   refD = getData("ref") 
   locnames <- NULL
   locstart <- 3 #row of locus start - 1
   tab3b[1,1] <- glabel(text="Profile:",container=tab3b)
   tab3b[2,1] <- glabel(text="Condition:",container=tab3b)
   tab3b[3,1] <- glabel(text="",container=tab3b)
   for(nm in 1:nM) { #print loci of each selected mixture
    tab3b[1,1 + nm] <- glabel(text=mixSel[nm],container=tab3b)
    tab3b[2,1 + nm] <- gcheckbox(text="",container=tab3b,checked=TRUE)
    subD <- mixD[[mixSel[nm]]] #select profile
    if(nm==1) locnames <- toupper(names(subD$adata)) #get loc-names
    for(i in 1:length(subD$adata)) { #use loci-name of first 
     locind <- grep(names(subD$adata)[i],locnames,fixed=TRUE) #get loc-index of stain:
     if(length(locind)==0) next #skip if not found
     if(nm==1) tab3b[locstart+i,1] <- locnames[i] #insert loc-name
     tab3b[locstart+locind,1 + nm]  <- gcheckbox(text="",container=tab3b,checked=TRUE)
     #note: some data may miss some loci. I.e. fine.
    }
   }
   if(nR>0) {
    for(nr in 1:nR) {
     tab3b[1,1 + nM + nr] <- glabel(text=refSel[nr],container=tab3b)
     tab3b[2,1 + nM + nr] <- gcheckbox(text="",container=tab3b,checked=FALSE)
     subD <- refD[[refSel[nr]]] #select profile
     for(i in 1:length(subD$adata)) {
      locind <- grep(names(subD$adata)[i],locnames,fixed=TRUE) #get loc-index of stain:
      if(length(locind)==0) next #skip if not found
      tab3b[locstart+locind,1 + nM + nr]  <- gcheckbox(text="",container=tab3b,checked=TRUE)
      #note: some data may miss some loci 
     }
    }
   } #end if refs.

  #User input info:
  tab3c[1,1] <- glabel(text="",container=tab3c)
  tab3c[2,1] <- glabel(text="#contributors:",container=tab3c)
  tab3c[2,2] <- gedit(text="2",container=tab3c,width=3)
  tab3c[3,1] <- glabel(text="Meth. param:",container=tab3c)
  tab3c[3,2] <- gedit(text="100",container=tab3c,width=3) #startvalue
  tab3c[4,1] <- glabel(text="Threshold:",container=tab3c)
  tab3c[4,2] <- gedit(text="50",container=tab3c,width=3)
  tab3c[5,1] <- glabel(text="Allow zeroMx:",container=tab3c)
  tab3c[5,2] <- gcheckbox(text="",container=tab3c,checked=FALSE)
 
  tab3c[9,1] <- glabel(text="",container=tab3c)  
  tab3c[10,1] = gbutton(text="Do deconvolution!",container=tab3c,handler=
   function(h,...) { 
    #take out selected data and send them to deconvolution
    mixD2 <- list() #will take data with locnames-order
    refD2 <- list() #will take data with locnames-order
    locsel_Mix <- numeric()
#Should be improved here: Latter stains should search correct loci.
    #1) Insert missing loci into ref-data:
    for(nm in 1:nM) { #for each selected mixes
     if(svalue(tab3b[2,1 + nm])) { #if checked
      mixD2[[mixSel[nm]]] = mixD[[mixSel[nm]]] #insert mixture
      locsel = rep(FALSE,length(locnames))
      for(i in 1:length(locnames)) { #for each needed locus
       if(svalue(tab3b[locstart+i,1 + nm])) locsel[i] = TRUE
      }
      locsel_Mix = cbind(locsel_Mix , locsel) #add column
     } #end if profile checked
    } #end for each mix
    if(length(mixD2)!=1) {
     tkmessageBox(message="Please select only one mix-profile! Multiple not implemented yet.")
     return
    }
 #Here we need to convert refData[[s]]$adata[[i]] to refData[[i]][[s]] for further use
    if(nR>0) { #note that some references may have other order of loci: This is taken into account.
     locsel_Ref <- matrix(FALSE,nrow=length(locnames),ncol=nR)
     for(i in 1:length(locnames)) { #for each needed locus
      refD2[[locnames[i]]] <- list()
      for(nr in 1:nR) { #for each selected references
       locind <- grep(locnames[i],toupper(names(refD[[refSel[nr]]]$adata)),fixed=TRUE) #get loc-index of prof
       locsel <- svalue(tab3b[locstart+locind,1 + nM + nr]) #boolean whether selected 
       if(length(locind)==0 | !locsel) {
        refD2[[locnames[i]]][[refSel[nr]]]  <- numeric() 
       } else {
         refD2[[locnames[i]]][[refSel[nr]]] <- refD[[refSel[nr]]]$adata[[locind]] #temp-store in list
         #only accept locus if checked
         if(svalue(tab3b[2,1 + nM + nr]) && locsel) locsel_Ref[locind,nr] = TRUE
       }
      } #end for each ref
     } #end for each locus     
    } #end if any refs

    #Call for Deconvolution: Rename ref-samples first:
     nCGUI <- as.numeric(svalue(tab3c[2,2])) #number of contributors
     rData <- NULL
     lsRef <- NULL
     condO <- NULL 
     if(nR>0 && any(locsel_Ref)) { #if any conditioned references
      rData<-refD2 
      lsRef<-locsel_Ref
      condO <- rep(0,nR)  #conditional order of the checked references
      pos = 1 #counter
      for(nr in 1:nR) { #for each selected references
       if(any(locsel_Ref[,nr])) {
        if(pos>nCGUI) tkmessageBox(message="Number of references exceeds number of contributors.")
        condO[nr] = pos #insert position
		pos = pos + 1
       }
      }
      print(locsel_Ref)
      print(rData)
     }
#only supports one stain!
      deconvlist <- deconvolve(nC=nCGUI,mixData=mixD2[[1]],refData=rData,condOrder=condO,locsel_Mix=locsel_Mix[,1],locsel_Ref=lsRef,eps=as.numeric(svalue(tab3c[3,2])),zeroMx=svalue(tab3c[5,2]),threshT=as.numeric(svalue(tab3c[4,2])),verbose=TRUE)
#     print(deconvlist)
     assign("deconvlist",deconvlist,envir=mmTK) #store results from deconv
     #send deconvolved results to next frame
     refreshTab4(layout="Layout2") #update result-tab
     svalue(nb) <- 4 #change notetab
   } #end handle function
  ) #end button
  tab3c[11,1] = glabel(text="",container=tab3c)  
 } #end refreshTab3

##############################################################
###############Tab 4: Deconvolution results:##################
##############################################################


 tab4tmp <- glayout(spacing=30,container=tab4)

 f_savetableDC = function(h,...) {
   deconvlist<-get("deconvlist",envir=mmTK) #load results from environment
   if(is.null(deconvlist)) {
    tkmessageBox(message="There is no deconvolution results available.")
    return
   }
   if(h$action=="Layout1") tab <- deconvlist$result1
   if(h$action=="Layout2") tab <- deconvlist$result2
   saveTable(addRownameTable(tab), "txt") #save 
  }

 refreshTab4 = function(layout="Layout2") {
   layouts <- c("Layout1","Layout2")
   deconvlist<-get("deconvlist",envir=mmTK) #load results from environment
   if(is.null(deconvlist)) return(1)
   if(layout=="Layout1") restab <- deconvlist$result1
   if(layout=="Layout2") restab <- deconvlist$result2
   tab4a = glayout(spacing=1,container=(tab4tmp[1,1] <-glayout(spacing=0,container=tab4tmp)),expand=TRUE) #table layout
   tab4b = glayout(spacing=1,container=(tab4tmp[2,1] <-glayout(spacing=0,container=tab4tmp)),expand=TRUE) #table of results
   tab4c = glayout(spacing=1,container=(tab4tmp[3,1] <-glayout(spacing=0,container=tab4tmp)),expand=TRUE) #storing result

   tab4b[1,1] <- gtable(addRownameTable(restab),container=tab4b,multiple=TRUE,width=mwW,height=mwH-2*mwH/3,do.autoscroll=TRUE,noRowsVisible=TRUE) #add to frame
   #create table in tab4a
   tab4b[1,1] <- glabel(text="",container=tab4b)
   tab4c[1,1] <- glabel(text="",container=tab4c)
   tab4c[2,1] <- glabel(text="             ",container=tab4c)
   tab4a[1,1] <- glabel(text="Table layout:",container=tab4a)
   tab4a[1,2] <- gradio(items=layouts,selected=which(layout==layouts),container=tab4a,horisontal=TRUE, handler=function(h,...) { refreshTab4( svalue(tab4a[1,2])  )})
   tab4c[2,2] <- gbutton(text="Save table",container=tab4c,handler=f_savetableDC,action=svalue(tab4a[1,2]))  

 }
 refreshTab4()

##############################################################
###############Tab 5: LR-calculation setup:###################
##############################################################

  #layout: 
  tab5tmp <- glayout(spacing=30,container=tab5)

   #function for assigning information from GUI to new:
   #overwrites existing information in "LRopt"
   getLRoptions = function(locnames,mixSel,refSel,tab5a,tab5b) { #send locnames,#mixes,#refs
    popLocs <- names(get("popFreq",envir=mmTK))  #get imported population frequencies. Needed for checking possible loci to calculate!
    nM <-length(mixSel)
    nR <-length(refSel)
    mixD = getData("mix")
    refD = getData("ref")
 
    LRopt <- list(mixData=list(),refData=list(),locnames=locnames) 
    for(nm in 1:nM) LRopt$mixData[[mixSel[nm]]] <- mixD[[mixSel[nm]]] #get selected data
    if(nR>0) {
     for(nr in 1:nR) LRopt$refData[[refSel[nr]]] <- refD[[refSel[nr]]] #get selected data
    }
    LRopt$Evidsel <- list()
    LRopt$Hp <- numeric()
    LRopt$Hd <- numeric()
    for(nm in 1:nM) { #get checked mixes and assign boolean for selected loci
     if(svalue(tab5a[1+nm,1])) LRopt$Evidsel[[mixSel[nm]]] <- rep(FALSE,length(mixD[[mixSel[nm]]]$adata))
    }
    if(nR>0) {
     for(nr in 1:nR) { #get checked refs
      if(svalue(tab5a[3+nM+nr,1])) LRopt$Hp <- c(LRopt$Hp,refSel[nr])
      if(svalue(tab5a[5+nM+nR+nr,1])) LRopt$Hd <- c(LRopt$Hd,refSel[nr])
     }
    }
    #get selected loci of checked mixtures
    for(nm in names(LRopt$Evidsel)) {  #for each selected evidence
     mixind <- grep(nm,mixSel,fixed=TRUE) #get mix-index in GUI
     subA <- names(mixD[[nm]]$adata) #get loci-names of selected mix-profile
     for(i in 1:length(subA)) { #for each locus in selected mix
      locind <- grep(subA[i],locnames,fixed=TRUE) #get loc-index of stain in GUI:
      if(length(locind)==0) next 
      if(length(grep(toupper(subA[i]),popLocs,fixed=TRUE))>0 && svalue(tab5b[1+locind,1 + mixind]) )  LRopt$Evidsel[[nm]][i] <- TRUE #locus was selected
     }
    }
    #get selected options:
    LRopt$uHp <- as.numeric(svalue(tab5a[9+nM+2*nR,2]))
    LRopt$uHd <- as.numeric(svalue(tab5a[10+nM+2*nR,2]))
    LRopt$DOprob <- as.numeric(svalue(tab5a[11+nM+2*nR,2]))
    LRopt$DIprob <- as.numeric(svalue(tab5a[12+nM+2*nR,2]))
    LRopt$theta <- as.numeric(svalue(tab5a[13+nM+2*nR,2]))
    LRopt$Qcalc <- svalue(tab5a[14+nM+2*nR,1]) #Qassignation
#    print(LRopt)
    return(LRopt)
   }

  

 #HELPFUNCTIONS CONSERNING LR calculations:
  calcLR <- function(LRopt,doLR=TRUE,verbose=TRUE) {
   #LRopt is first a object from LRoption-function
   #if doLR is FALSE, the LRopt-object is returned without LR-calculations
   require(forensim) #this calculation requires the forensim package
#   LRopt<- get("LRopt",envir=mmTK)  #get imported population frequencies
   popFreq <- get("popFreq",envir=mmTK)  #get imported population frequencies
   minFreq <- getminFreq(get("freqSize",envir=mmTK)) #get assigned minFrequency
   LRfit <- list() #store organized data

   #organize data:
   maxC <- max(LRopt$uHp+length(LRopt$Hp),LRopt$uHd+length(LRopt$Hd))

   #if dropout-is specified in action-list:
   if(length(LRopt$DOprob)>1 | length(LRopt$DIprob)>1) print("Calculating LR for different prD- and prC-values")

   #Find loci in population frequencies: 
   #NB: Data-loci must be a subset of popfreq-loci!
   #If some loci missing in one of the reference: The loci is not calculated!
   evidnames <- names(LRopt$Evidsel)
   nM <- length(evidnames)
   for(ln in LRopt$locnames) { #for each locus
    poplocind <-  grep(toupper(ln), toupper(names(popFreq)),fixed=TRUE)
    if(length(poplocind)==0) next #skip not-existing in pop
    if(verbose) print(paste("Calculating LR for loci ",ln,"...",sep=""))
    freq <- popFreq[[poplocind]] #take out frequeny
    if(is.null(freq)) next #no frequencies found for given allele
    freqQ <- freq #default it is all alleles
    LRfit[[ln]] <- list() #storage for each loci
    evidA <- numeric() #alleles for evidence
    evidList <- list()
    for(nm in evidnames) { #for each evidence
     evidlocind <- grep(ln,names(LRopt$mixData[[nm]]$adata),fixed=TRUE) #get lociind in evidence
     if(length(evidlocind)==0 || !LRopt$Evidsel[[nm]][evidlocind]) {
      tmpA=0 #locOK <- FALSE #loci was not found in evidence or not selected
      evidList[[nm]] <- numeric()
     } else {
      tmpA <- LRopt$mixData[[nm]]$adata[[evidlocind]]
      if(length(tmpA)==0) tmpA <- 0 #set to zero of no contributors
      evidList[[nm]] <- LRopt$mixData[[nm]]$adata[[evidlocind]]
     } 
     if(length(evidA)==0) { evidA <- tmpA #add new
     } else { evidA <- c(evidA,0,tmpA) } #add to existing: REQUIRED IN likEvid()
    }
  
    #get checked references 
    refHp <- as.numeric()
    nonContHd <- as.numeric() #references in Hp but not in Hd: sot
    hpnames <- as.numeric()
    nchdnames <- as.numeric()
    for(nr in LRopt$Hp) { #for each reference in Hp
     reflocind <- grep(ln,names(LRopt$refData[[nr]]$adata),fixed=TRUE) #get loc-index of reference:
     if(length(reflocind)==0 || length(LRopt$refData[[nr]]$adata[[ln]])==0) {
      print(paste("Hp: Loci",ln,"was not found in reference ",nr,". You should unselect loci."))      
     } else {
      refHp <- cbind(refHp , LRopt$refData[[nr]]$adata[[ln]]) #get reference-data
      hpnames <- c(hpnames,nr) 
      if(!any(nr==LRopt$Hd)) {
       nonContHd <- cbind(nonContHd , LRopt$refData[[nr]]$adata[[ln]]) 
       nchdnames <- c(nchdnames,nr)
      }
     }
    }
    if(length(hpnames)>0) colnames(refHp) <- hpnames
    if(length(nchdnames)>0) colnames(nonContHd) <- nchdnames
    hdnames <- as.numeric()
    refHd <- as.numeric()
    for(nr in LRopt$Hd) { #for each reference in Hp
     reflocind <- grep(ln,names(LRopt$refData[[nr]]$adata),fixed=TRUE) #get loc-index of reference:
     if(length(reflocind)==0 || length(LRopt$refData[[nr]]$adata[[ln]])==0) {
      print(paste("Hd: Loci",ln,"was not found in reference ",nr,". You should unselect loci."))      
     } else {
      refHd <- cbind(refHd , LRopt$refData[[nr]]$adata[[ln]]) #get reference-data
      hdnames <- c(hdnames,nr) 
     }
    }  #END of data organization
    if(length(hdnames)>0) colnames(refHd) <- hdnames

#print(ln)
#print(evidA)
#print(refHp)
#print(refHd)

    #frequency-handling:
    allA <- unique(c(evidA,refHp,refHd)) 
    allA <- allA[allA!=0] #not zeros
    freqN <- names(freq)
    #insert missing allele in 'freq' if some are missing:
    newA <- allA[!allA%in%freqN] #new alleles
    if(length(newA)>0) {
     print(paste0("In locus ",ln,": Allele(s) ",newA," was inserted with min. frequency ",prettyNum(minFreq)))
     freq <- c(freq,rep(minFreq,length(newA)))
     freq <- freq/sum(freq) #normalize
     names(freq) <- c(freqN,newA)
     popFreq[[poplocind]] <- freq #update frequency table
    }
    freqQ  <- freq #default is all frequencies

    #if Q-assignated: 
    if(LRopt$Qcalc) {
      freqQ <- freq[freqN%in%allA] #truncate alleles
      freqQ <- c(freqQ,1-sum(freqQ))
      names(freqQ)[length(freqQ)] = "99"
    }
    #start calculate LR:
    nPrD <- length(LRopt$DOprob)
    nPrC <- length(LRopt$DIprob)
    numM <- denoM <- lrM <- matrix(1,nrow=nPrD,ncol=nPrC)
    colnames(numM) <- colnames(denoM) <- colnames(lrM) <- LRopt$DIprob
    rownames(numM) <- rownames(denoM) <- rownames(lrM) <- LRopt$DOprob
    for(ddin in 1:nPrC) { #for each dropin
     for(ddout in 1:nPrD) { #for each dropout
      if(verbose && (ddout%%5==0)) { #print for each 5. calculation
       print(paste( ((ddin-1)*nPrD + ddout - 1) / (nPrD*nPrC)*100, "% complete",sep=""))
      }
      PrD <- rep(LRopt$DOprob[ddout],maxC)
      PrC <- rep(LRopt$DIprob[ddin],maxC) 
      numM[ddout,ddin] <- denoM[ddout,ddin] <- lrM[ddout,ddin] <- 1  #default value
      if(length(evidA)>0 && !all(evidA==0) && doLR) { #if any samples to calculate
       lrobj <- LR(evidA,Tp=c(refHp),Td=c(refHd),Vp=NULL,Vd=c(nonContHd), xp=LRopt$uHp, xd=LRopt$uHd, theta=LRopt$theta, prDHet=PrD, prDHom=PrD^2, prC=PrC, freq=freqQ)  
       numM[ddout,ddin] <- lrobj$num
       denoM[ddout,ddin] <- lrobj$deno
       lrM[ddout,ddin] <- lrobj$LR
      }
      #store to LRfit-object: Loci-results given
      LRfit[[ln]] <- list(hp=numM,hd=denoM,LR=lrM,evidList=evidList,evid=evidA,refHp=refHp,refHd=refHd,nonContHd=nonContHd,freq=freq,freqQ=freqQ)
     } #end for each dropin
    } #end for each dropout
   } #end for each loci

   assign("popFreq",popFreq,envir=mmTK)  #assigning updated popFreq
   if(!doLR) return(LRfit)
   LRopt$LRfit <- LRfit
   assign("LRopt",LRopt,envir=mmTK) #store selected LR-model
   if(nPrC==1 & nPrD==1) { #if ordinary LR-calculation
    refreshTab6() #refresh to result-table 
   } else { #if precalculation
    #plot LR as a function of prD and prC:
    locs <- names(LRfit)
    jointlr = 1
    for(loc in locs) {
     jointlr <- jointlr*LRfit[[loc]]$LR
    }
    if(ncol(jointlr)==1) { #if only drop-out calc:
     pdx <- as.numeric(rownames(jointlr))
     plot(pdx,log10(jointlr),xlab="Pr_D",ylab="log_10(LR(D))",main="Sensitivity plot")
     lines(spline(pdx,log10(jointlr),n=100))
    } else {
     print("multiple drop-in probabilities not implemented yet!")
    }
   } #end if precalc
  } #end LR-calculation function
 ####### END CALC LR#######

   
   #precalcLR function of dropout
   precalcLR = function(LRopt,tab5c) { 
    xx <- as.numeric(unlist(strsplit(svalue(tab5c[4,2]),"-")))
    prDvec <- seq(xx[1],xx[2],l=as.numeric(svalue(tab5c[5,2])))
    LRopt$DOprob <- prDvec #change drop-out range
    calcLR(LRopt,doLR=TRUE,verbose=TRUE) #call LR-function
   }

   precalcDO_CI = function(LRopt,tab5c,type="weight") {
     #type="weight": calculate DO-distr under both hypothesis
     #type="dbsearch": calculate DO-distr under hd only
     LRopt$LRfit <- calcLR(LRopt,doLR=FALSE,verbose=FALSE) #get need info of LRcalculation 
     LRopt$dropquant <- list()
     sn <- names(LRopt$Evidsel)
     alph <- as.numeric(svalue(tab5c[8,2]))
     qqs <- c(alph,0.5,1-alph)
     par(mfrow=c(length(sn),1))
     for(ss in sn) {
      droplist <- calcDOdistr(LRopt,sample=ss,minS=as.numeric(svalue(tab5c[9,2])),hdOnly=(type=="dbsearch")) #sample only hd-case if dbsearch
      plot(density(droplist$hdsamples,from=0,to=1),main=paste("DO-distr for sample",ss))
      hpCI <- ""
      if(type=="weight") {
       qqHp <- round(quantile(droplist$hpsamples,qqs),3)
       print("DO-quantiles for Hp:") 
       print(qqHp)
       LRopt$droplist[[ss]] <- qqHp
       lines(density(droplist$hpsamples,from=0,to=1),lty=2)
       for(kk in 1:length(qqHp)) lines(c(qqHp[kk],qqHp[kk]),c(0,100),lty=2,lwd=0.7,col=2)
       hpCI <- paste0("Hp[",alph/2,",0.5,",1-alph/2,"]=[",qqHp[1],",",qqHp[2],",",qqHp[3],"]")
      }
      qqHd <- round(quantile(droplist$hdsamples,qqs),3)
      print("DO-quantiles for Hd:") 
      print(qqHd)
      LRopt$droplist[[ss]] <- qqHd
      for(kk in 1:length(qqHd)) lines(c(qqHd[kk],qqHd[kk]),c(0,100),lwd=0.7,col=2)
      hdCI <- paste0("Hd[",alph/2,",0.5,",1-alph/2,"]=[",qqHd[1],",",qqHd[2],",",qqHd[3],"]")
      legend("topright",c(hdCI,hpCI),lty=1:2)
     }
     par(mfrow=c(1,1))
   }

  refreshTab5 = function(mixSel,refSel,type="weight") { 
   #type={"weight","dbsearch"}
   #weight: must have both mixture and reference profiles
   #dbsearch: must have both mixture and database, reference profiles is optional
   #IF LRopt=NULL; get assigned values?
   popLocs <- names(get("popFreq",envir=mmTK))  #get imported population frequencies. Needed for checking possible loci to calculate!
   mixD = getData("mix")
   refD = getData("ref") 
   nM = length(mixSel) #number of mix-profiles
   nR = length(refSel) #number of ref-profiles
   locnames <- NULL
   locstart <- 3 #row of locus start - 1

   tab5a = glayout(spacing=0,container=(tab5tmp[1,1] <-gframe("Model configuration",container=tab5tmp))) 
   tab5b = glayout(spacing=0,container=(tab5tmp[1,2] <-gframe("Data selection",container=tab5tmp))) 
   tab5c = glayout(spacing=0,container=(tab5tmp[1,3] <-gframe("Precalculations",container=tab5tmp))) 
   tab5d = glayout(spacing=0,container=(tab5tmp[2,2] <-gframe("Calculations",container=tab5tmp)))  

   #Hypothesis selection:
   tab5a[1,1] <- glabel(text="Evidence(s):",container=tab5a)
   for(nm in 1:nM) {
    tab5a[1+nm,1] <- gcheckbox(text=mixSel[nm],container=tab5a,checked=TRUE)
   }
   tab5a[2+nM,1] <- glabel(text="",container=tab5a) #space
   if(type=="weight") tab5a[3+nM,1] <- glabel(text="Contributor(s) under Hp:",container=tab5a)
   if(type=="dbsearch") tab5a[3+nM,1] <- glabel(text="Contributor(s) under Hp\n (DB-reference already included):",container=tab5a)
   tab5a[4+nM+nR,1] <- glabel(text="",container=tab5a) #space
   tab5a[5+nM+nR,1] <- glabel(text="Contributor(s) under Hd:",container=tab5a)
   tab5a[6+nM+2*nR,1] <- glabel(text="",container=tab5a) #space
   tab5a[7+nM+2*nR,1] <- glabel(text="",container=tab5a) #space

   if(nR>0) { #if any references
    for(nr in 1:nR) {
     tab5a[3+nM+nr,1] <- gcheckbox(text=refSel[nr],container=tab5a,checked=TRUE)
     tab5a[5+nM+nR+nr,1] <- gcheckbox(text=refSel[nr],container=tab5a,checked=TRUE)
    }
   }
 
   #Parameter selection:
   tab5a[8+nM+2*nR,1] <- glabel(text="Parameters:",container=tab5a)
   tab5a[9+nM+2*nR,1] <- glabel(text="unknowns (Hp): ",container=tab5a)
   tab5a[9+nM+2*nR,2] <- gedit(text="1",container=tab5a,width=4)
   tab5a[10+nM+2*nR,1] <- glabel(text="unknowns (Hd): ",container=tab5a)
   tab5a[10+nM+2*nR,2] <- gedit(text="2",container=tab5a,width=4)
   tab5a[11+nM+2*nR,1] <- glabel(text="Probability of Drop-out: ",container=tab5a)
   tab5a[11+nM+2*nR,2] <- gedit(text="0.1",container=tab5a,width=4)
   tab5a[12+nM+2*nR,1] <- glabel(text="Probability of Drop-in: ",container=tab5a)
   tab5a[12+nM+2*nR,2] <- gedit(text="0.05",container=tab5a,width=4)
   tab5a[13+nM+2*nR,1] <- glabel(text="Theta-correction: ",container=tab5a)
   tab5a[13+nM+2*nR,2] <- gedit(text="0",container=tab5a,width=4)
   tab5a[14+nM+2*nR,1] <- gcheckbox(text="Q-assignation",container=tab5a,checked=FALSE,horisontal=TRUE)

   #note: User includes those who are interest to select for conditioning directly!
   #Locus selecter (same as in deconvolution)
   #'locinames' - assigned as first mixture data
   #Locinames in References are then matched with 'locinames'
   #user may select loci of mixtures 
   tab5b[1,1] <- glabel(text="Loci:",container=tab5b)
   for(nm in 1:nM) {
    tab5b[1,1 + nm] <- glabel(text=mixSel[nm],container=tab5b)
    subD <- mixD[[mixSel[nm]]] #select profile
    if(nm==1) locnames <- toupper(names(subD$adata)) #get loc-names from first mixture!
    for(i in 1:length(subD$adata)) { #for each locus
     locind <- grep(names(subD$adata)[i],locnames,fixed=TRUE) #get loc-index of stain:
     if(nm==1) tab5b[1+i,1] <- locnames[i] #insert loc-name
     tab5b[1+locind,1 + nm]  <- gcheckbox(text="",container=tab5b,checked=TRUE)
     if(length(grep(names(subD$adata)[i],popLocs,fixed=TRUE))==0) enabled(tab5b[1+locind,1 + nm]) <- FALSE #deactivate non-existing locus
    }
   }  

   if(nR>0) { #if any references 
    #show referenced (possible partial profiles) profiles
    for(nr in 1:nR) { #for each reference
     tab5b[1,1 + nM + nr] <- glabel(text=refSel[nr],container=tab5b) #name of reference
     subD <- refD[[refSel[nr]]] #select profile
     for(i in 1:length(subD$adata)) { #for each marker 
      locind <- grep(toupper(names(subD$adata))[i],locnames,fixed=TRUE) #get loc-index of stain:
      if(length(locind)==0) next #skip marker if not existing
      check <- TRUE
      if(length(subD$adata[[i]])==0) check <- FALSE
      tab5b[1+locind,1 + nM + nr]  <- gcheckbox(text="",container=tab5b,checked=check)
      enabled(tab5b[1+locind,1 + nM + nr]) <- FALSE #they cannot be selected!!
      #note: some data may miss some loci 
     }
    }
   }

   #add buttons:
   if(type=="weight") {
    tab5c[1,1] = glabel(text="",container=tab5c)
    tab5c[2,1] = glabel(text="",container=tab5c)
    tab5c[3,1] = glabel(text="Drop-out plot: ",container=tab5c)
    tab5c[3,2] = gbutton(text="Calc!",container=tab5c,handler=
     function(h,...) {
     #send GUI-objects to get variables in LRopt which is sent to precalcLR together with GUI option
     precalcLR(getLRoptions(locnames,mixSel,refSel,tab5a,tab5b),tab5c) #
    })   
    tab5c[4,1] = glabel(text="Range: ",container=tab5c)
    tab5c[4,2] = gedit(text="0.1-0.6",container=tab5c,width=8) #range of dropout
    tab5c[5,1] = glabel(text="#ties=",container=tab5c)
    tab5c[5,2] = gedit(text="9",container=tab5c,width=3) #number of ties
   }
   tab5c[6,1] = glabel(text="",container=tab5c)
   tab5c[7,1] = glabel(text="Drop-out distr:",container=tab5c)
   tab5c[7,2] = gbutton(text="Calc!",container=tab5c,handler= 
    function(h,...) {
     #send GUI-objects to get variables in LRopt which is sent to calcDO_CI together with GUI option
     precalcDO_CI(getLRoptions(locnames,mixSel,refSel,tab5a,tab5b),tab5c,type=type)
    })
   tab5c[8,1] = glabel(text="alpha=",container=tab5c)
   tab5c[8,2] = gedit(text="0.05",container=tab5c,width=4) #alpha
   tab5c[9,1] = glabel(text="#Samples=",container=tab5c)
   tab5c[9,2] = gedit(text="500",container=tab5c,width=4) #number of samples

  if(type=="weight") {
    tab5d[1,1] = gbutton(text="Do LR calculation!",container=tab5d,handler=
	function(h,...) {
     #send GUI-objects to get variables in LRopt which is sent to calcLR
     calcLR(getLRoptions(locnames,mixSel,refSel,tab5a,tab5b),verbose=FALSE)  
     svalue(nb) <- 6 #change notetab
   })
  }
  if(type=="dbsearch") {
    tab5d[1,1] = gbutton(text="Do database search!",container=tab5d,handler=
	function(h,...) {
     #send GUI-objects to get variables in LRopt which is sent to calcLR     
     doDBsearch(getLRoptions(locnames,mixSel,refSel,tab5a,tab5b),verbose=FALSE)  
     refreshTab7() #update result-tab
     svalue(nb) <- 7 #change notetab
   })
  }


  } #end refresh table 

  #options:
  #Q-assignation
  #number of contr


#############################################################
#############Tab 6: LR-calculation results:##################
#############################################################

  tab6tmp <- glayout(spacing=30,container=tab6)

 f_savetableLR = function(h,...) {
   LRopt<-get("LRopt",envir=mmTK) #load results from environment
   if(is.null(LRopt$LRtab) || is.null(LRopt$joint)) {
    tkmessageBox(message="There is no LR-results available.")
    return
   }
   tab <- LRopt$LRtab
   tab <- rbind(tab,LRopt$joint)
   saveTable(addRownameTable(tab),"txt")
 }

  refreshTab6 = function(acc=5) { 
   LRopt<-get("LRopt",envir=mmTK) #load weight-of-evidence results from environment
   if(is.null(LRopt) || is.null(LRopt$LRfit)) return(1);
   #acc is number of rounding
   tab6a = glayout(spacing=0,container=(tab6tmp[1,1] <-glayout(spacing=0,container=tab6tmp))) 
   tab6b = glayout(spacing=0,container=(tab6tmp[1,2] <-glayout(spacing=0,container=tab6tmp))) 
   tab6c = glayout(spacing=0,container=(tab6tmp[2,1] <-glayout(spacing=0,container=tab6tmp))) 
   tab6d = glayout(spacing=0,container=(tab6tmp[2,2] <-glayout(spacing=0,container=tab6tmp))) 
   LRfit <- LRopt$LRfit
   nL <- length(LRfit) #number of loci:
   locnames <- as.numeric()
   hpJoint <- as.numeric()
   hdJoint <- as.numeric()
   for(ln in names(LRfit)) {
    locnames <- c( locnames,ln)
    hpJoint <- c(hpJoint,LRfit[[ln]]$hp )
    hdJoint <- c(hdJoint,LRfit[[ln]]$hd )
   }
   lrJoint<- hpJoint/hdJoint
   LRopt$LRtab <- cbind(hpJoint,hdJoint,lrJoint)
   rownames(LRopt$LRtab) <- locnames
   LRopt$joint <- c(prod(hpJoint),prod(hdJoint),prod(lrJoint))
   LRopt$joint <- rbind(LRopt$joint,log10(LRopt$joint))
   colnames(LRopt$joint) <- colnames(LRopt$joint) <- c("Hp","Hd","LR")
   rownames(LRopt$joint) <- c("Joint","log10Joint")
   assign("LRopt",LRopt,envir=mmTK) #store/update results

   rr <- function(x) signif(x,acc) #user may specify
   LRtab <-  cbind(rownames(LRopt$LRtab),rr(LRopt$LRtab)) #format table for outprint
   colnames(LRtab) <- c("loci","Hp","Hd","LR")
   tab6b[1:nrow(LRtab),1:4] <- gtable(LRtab,container=tab6b,multiple=FALSE,height=min(nL*22,500)) #plot LR-results
   tab6b[nrow(LRtab)+1,2] <- glabel("",container=tab6b)
   tab6b[nrow(LRtab)+2,2] <- glabel("Hp",container=tab6b)
   tab6b[nrow(LRtab)+2,3] <- glabel("Hd",container=tab6b)
   tab6b[nrow(LRtab)+2,4] <- glabel("LR",container=tab6b)
   tab6b[nrow(LRtab)+3,1] <- glabel("Joint",container=tab6b)
   tab6b[nrow(LRtab)+4,1] <- glabel("log10Joint",container=tab6b)
  
   tab6b[nrow(LRtab)+3,2] <- glabel(rr(LRopt$joint[1,1]),container=tab6b)
   tab6b[nrow(LRtab)+3,3] <- glabel(rr(LRopt$joint[1,2]),container=tab6b)
   tab6b[nrow(LRtab)+3,4] <- glabel(rr(LRopt$joint[1,3]),container=tab6b)
   tab6b[nrow(LRtab)+4,2] <- glabel(rr(LRopt$joint[2,1]),container=tab6b)
   tab6b[nrow(LRtab)+4,3] <- glabel(rr(LRopt$joint[2,2]),container=tab6b)
   tab6b[nrow(LRtab)+4,4] <- glabel(rr(LRopt$joint[2,3]),container=tab6b)

   printModel(LRopt,guitab=tab6a) #print LR-model

   #options:
   tab6c[1,1] <- glabel("Tippet:",container=tab6c)
   tab6c[2,1] <- gbutton(text="RM precalc",container=tab6c,handler=
    function(x) {
     LRopt$LRRMlist <- getLRRMlist(LRopt,svalue(tab6a[4+length(LRopt$Evidsel),1]))
     assign("LRopt",LRopt,envir=mmTK) #store LRRM-info
     refreshTab6(acc=as.numeric(svalue(tab6d[2,2]))) #refresh tab
    })
   tab6c[2,2] <- gbutton(text="Pval RM",container=tab6c,handler=
    function(x) {
      print(LRopt$joint[1,3]) #observed LR
      pval <- calcPvalue(LRopt$LRRMlist,LRopt$joint[1,3]) #send random man calculations and observed LR
      msg <- paste("Calculated pvalue: ",pval,sep="")
      print(msg)
      gmessage(message=msg,title="Calculate p-value",icon="info")
    }
   )
   tab6c[3,2] <- gedit("1e6",container=tab6c,width=4)
   tab6c[3,1] <- gbutton(text="Tippet plot",container=tab6c,handler=
     function(x) {
      lr0 <- LRopt$joint[2,3] #get observed log10-LR
      M  <- as.numeric(svalue(tab6c[3,2]))
      RM_LR = 1
      Glist = LRopt$LRRMlist$Glist #stored genotype-information
      print("Sampling tippets...")
      for(i in 1:length(Glist)) {
       if(is.null(Glist[[i]])) next
       X = sample(1:length(Glist[[i]]$Gprob),M,replace=TRUE,prob=Glist[[i]]$Gprob)
       RM_LR = RM_LR*LRopt$LRRMlist$LRRM[[i]][X]
      }
      RM_LR = log10(RM_LR)
      minmax = range(RM_LR)
      qvals <- c(0.01,0.05,0.5,0.95,0.99)
      quantiles<-quantile(RM_LR,qvals)
      tab <- cbind(c("min",as.character(qvals),"max"),round(c(minmax[1],quantiles,minmax[2]),4))
      colnames(tab) <- c("qq","log10LR(qq)")
      print(tab)
      print(paste0("Discrimanatory metric(99% quantile)=",round(lr0-quantiles[length(qvals)],4))) 
      #show table in window:
      tipquantwin <- gwindow("Tippet quantiles")
      gtable(as.data.frame(tab),container=tipquantwin,width=200,height=170)

      print("...finished. Now plotting tippet-distribution.")
      mtxt = paste("Tippet calculation with M=",M," iterations.",sep="")
      plot(ecdf(RM_LR),xlim=c(minmax[1],max(minmax[2],lr0)) , main=mtxt,xlab="log10(LR)")
      points(lr0,1,pch=10,col="blue")
      lines(rep(lr0,2),c(1,0),lty=1,col="blue",lwd=0.5)

     })
   if(is.null(LRopt$LRRM)) enabled(tab6c[2,2]) <- FALSE
   if(is.null(LRopt$LRRM)) enabled(tab6c[3,1]) <- FALSE

   tab6d[1,1] <- glabel("Options:",container=tab6d)
   tab6d[2,1] <- glabel("Rounding decimal:",container=tab6d)
   tab6d[2,2] <- gedit(paste(acc),container=tab6d,width=3)
   tab6d[3,1] <- gbutton(text="Refresh table",container=tab6d,handler=
    function(h,...) {
    refreshTab6(acc=as.numeric(svalue(tab6d[2,2])))
   } )  
   tab6d[3,2] <- gbutton(text="Save table",container=tab6d,handler=f_savetableLR)  
  }
  refreshTab6()

##############################################################
###############Tab 7: Database search:########################
##############################################################

  tab7tmp <- glayout(spacing=30,container=tab7) #grid-layout

 f_savetableDS = function(h,...) {
   DBsearch <-get("DBsearch",envir=mmTK) #load results from environment
   if(is.null(DBsearch)) {
    tkmessageBox(message="There is no database search results available.")
    return
   }
   sep="txt"
   tabfile = gfile(text="Save table",type="save")
   if(!is.na(tabfile)) {
    if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
    if(h$action=="LR") ord <- order(as.numeric(DBsearch[,2]),decreasing=TRUE) 
    if(h$action=="MAC") ord <- order(as.numeric(DBsearch[,3]),decreasing=TRUE) 
    write.table(DBsearch[ord,],file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) #load environment
    print(paste("Full ranked table saved in ",tabfile,sep=""))
   }
  }

 refreshTab7 = function(layout="LR") {
   DBsearch <- get("DBsearch",envir=mmTK) #load results from environment
   LRoptDB <- get("LRoptDB",envir=mmTK)
   if(is.null(DBsearch) || is.null(LRoptDB) ) return(1);
   tab7a = glayout(spacing=1,container=(tab7tmp[1,1] <-glayout(spacing=0,container=tab7tmp)),expand=TRUE) #table layout (sorted)
   tab7b = glayout(spacing=1,container=(tab7tmp[1,2] <-glayout(spacing=0,container=tab7tmp)),expand=TRUE) #storing result
   tab7c = glayout(spacing=1,container=(tab7tmp[2,1] <-glayout(spacing=0,container=tab7tmp)),expand=TRUE) #table of model
   tab7d = glayout(spacing=1,container=(tab7tmp[2,2] <-glayout(spacing=0,container=tab7tmp)),expand=TRUE) #table of results

   layouts <- c("LR","MAC")
   tab7a[1,1] <- glabel(text="Table sort:",container=tab7a)
   tab7a[1,2] <- gradio(items=layouts,selected=grep(layout,layouts,fixed=TRUE),container=tab7a,horisontal=TRUE, handler=function(h,...) { refreshTab7( svalue(tab7a[1,2])  )})
   tab7b[1,1] <- glabel(text="",container=tab7c)
   tab7b[2,1] <- glabel(text="             ",container=tab7b)
   tab7b[2,2] <- gbutton(text="Save table",container=tab7b,handler=f_savetableDS,action=svalue(tab7a[1,2]))  
   printModel(LRoptDB,guitab=tab7c) #print LR-model

   if(layout=="LR") ord <- order(as.numeric(DBsearch[,2]),decreasing=TRUE) 
   if(layout=="MAC") ord <- order(as.numeric(DBsearch[,3]),decreasing=TRUE) 
   if(length(ord)>10000) {
    print("Database contained more than 10000 samples.\n Showing first 10000 ranked samples only!")
    ord <- ord[1:10000] #only first 10000 shown
   }
   tab7d[1,1] <- gtable(DBsearch[ord,] ,container=tab7d,multiple=TRUE,height=mwH-2*mwH/3,do.autoscroll=TRUE,noRowsVisible=TRUE) #add to frame
 }
 refreshTab7()

 visible(mainwin) <- TRUE


} #end funcions
