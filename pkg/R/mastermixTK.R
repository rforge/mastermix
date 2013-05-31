###########
#Changelog#
###########
#14.05.13 - Start create GUI for mastermix-function

#' @title mastermixTK
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description mastermixTK is a GUI for the function mastermix for performing linear deconvolution procedure of STR DNA mixtures.
#' @export
#' @usage mastermixTK()
#' @details The function is a graphical layer for the function mastermix. See ?mastermix for more information.
#' 
#' The structure of profiles is: Data[['samplename']]$adata[['locusname']] for allele information and Data[['samplename']]$hdata[['locusname']] for allele height information (if any)
#' 
#' @keywords deconvolution, optimization

#rm(list=ls()) #must be removed after
# require(mastermix)

mastermixTK = function() {
 source( system.file("getKit_05.R",package="mastermix") )
 source( system.file("plotEPG.R",package="mastermix") )

 options(guiToolkit="tcltk")

 #Required in GUI:
 require(gWidgetstcltk) #requires only gWidgets also:

 #Required in mastermix-function:
 require(gtools)
 require(MASS)

 #version:
 version = 1

 #create own environment within mastermixTK - function (i.e. env-object is enclosed from this function)
 mmTK = new.env()
 #assign("name",1,envir=mmTK) #assign to environment
 #get(name,envir=mmTK) #grab from environment
 #save(mmTK,filename="projectname.Rdata") #can save environment to file

 #initialize environment variables
 assign("workdir",NULL,envir=mmTK) #assign working directory to mmTK-environment
 assign("freqfolder",NULL,envir=mmTK) #assign freqfolder to mmTK-environment
 assign("mixData",NULL,envir=mmTK) #assign mixdata to mmTK-environment
 assign("refData",NULL,envir=mmTK) #assign refdata to mmTK-environment
 assign("kits",NULL,envir=mmTK) #assign kitList to mmTK-environment
 assign("popFreq",NULL,envir=mmTK) #assign popFreq to mmTK-environment
 assign("deconvlist",NULL,envir=mmTK) #assign popFreq to mmTK-environment

 #Function to get data from environment
 getData = function(type) {
   Data <- NULL
   if(type=="mix") Data <- get("mixData",envir=mmTK) #assign kit to mmTK-environment
   if(type=="ref") Data <- get("refData",envir=mmTK) #assign kit to mmTK-environment 
   return(Data)
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
   for(j in 1:ncol(tab)) {
     tmp = tab[,j]
     tmp2 = tmp[!is.na(tmp)]
     names(tmp2) = Anames[!is.na(tmp)]
     freqlist[[j]] = tmp2
   }
   names(freqlist) = colnames(tab)
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

 #Robust function for reading tables:
 tableReader=function(filename) {
  tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE)
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE)
  return(tab)
 }


 strsplit2 <- function(x,spl) {
  if(nchar(x)==0) return("")
  txt <- x
  for(j in 1:length(spl)) {
   txt <- unlist(strsplit(txt,split=spl[j]))
  }
  return(txt)
 }
 
 #########################
 ######DATA LOADING#######
 #########################


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
  lind = grep("marker",tolower(cn)) #locus col-ind
  sind = grep("sample",tolower(cn)) #sample col-ind
  A_ind = grep("allele",tolower(cn)) #allele col-ind
  H_ind = grep("height",tolower(cn)) #height col-ind
  ln = unique(X[,lind]) #locus names
  sn = unique(X[,sind]) #sample names
  I = length(ln)
  Y = list() #insert non-empty characters:
  for(k in 1:length(sn)) { #for each sample
   Y[[sn[k]]] = list() #one list for each sample
   if(length(A_ind)>0) Y[[sn[k]]]$adata=list()
   if(length(H_ind)>0) Y[[sn[k]]]$hdata=list()
   for(i in 1:I) { #for each locus
     xind = X[,sind]==sn[k] & X[,lind]==ln[i] #get index in X for given sample and locus
     if(length(A_ind)>0) Y[[sn[k]]]$adata[[ln[i]]] = as.character(X[xind,A_ind][!X[xind,A_ind]%in%c("","NA")])
     if(length(H_ind)>0) {
      Y[[sn[k]]]$hdata[[ln[i]]] = as.numeric(as.character(X[xind,H_ind][!X[xind,H_ind]%in%c("","NA")]))
    }
   }
  }
  return(Y)
 }

 ########################
 #######VIEW DATA########
 ########################

#wrapper function for plotting mixture data as epg profile
plotEPG <- function(Data,kit,sname="") {
 #mixData is list with allele and height data. Only one sample!
 #for selected sample:
 if(all(length(unlist(Data$hdata))==0)) { 
  tkmessageBox(message="There is no height data in sample.")
 } else if(all(length(unlist(Data$adata))==0)) {
  tkmessageBox(message="There is no allele data in sample!")
 } else {
  generateEPG(typingKit=kit,alleleList=Data$adata,peakHeightList=Data$hdata, locusVector=names(Data$adata),sampleName=sname, drawBoxPlots=FALSE, drawPeaks=TRUE)
 }
}

###################################################################
###########################GUI#####################################
###################################################################

###########GUI WINDOW STARTS#########

 #Menu bar file-lists:
 f_setwd = function(h,...) {
  dirfile = gfile(text="Select folder",type="selectdir")
  if(!is.na(dirfile)) {
   setwd(dirfile)
   assign("workdir",dirfile,envir=mmTK) #assign working directory
  }
 }
 f_openproj = function(h,...) {
  projfile = gfile(text="Open project",type="open")
  if(!is.na(projfile)) load(projfile) #load environment
 }
 f_saveproj = function(h,...) {
  projfile = gfile(text="Save project",type="save")
  if(!is.na(projfile)) {
   save(mmTK,file=projfile) #save environment
   print(paste("Project saved in ",projfile,sep=""))
  }
 }
 f_quitproj = function(h,...) {
  ubool <- gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(svalue(ubool)) {
   savefile <- gfile(text="Save file",type="save")
   save(mmTK,file=savefile) 
   print(paste("Project saved in ",savefile,sep=""))
  } else { 
   print("Program terminated without saving")
  }
  dispose(mainwin) #remove window
 }
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
 
 #tmp
 #User should be able to select kit and population in tab
#  kitSel=1
#  popSel=1
 

 #Main window:
 mainwin <- gwindow(paste("Mixture Interpretation Tool v",version,sep=""), visible=FALSE)
 gmenu(mblst,container=mainwin)
 nb = gnotebook(container=mainwin)
 tab1 = glayout(spacing=100,container=nb,label="Create data")
 tab2 = glayout(spacing=100,container=nb,label="Import data")
 tab3 = glayout(spacing=100,container=nb,label="Deconvolution Setup")
 tab4 = glayout(spacing=100,container=nb,label="Deconvolution Results",expand=TRUE) 
 svalue(nb) <- 2 #initial start at second tab

####################################################
###############Tab 1: Create Data:##################
####################################################
 #layout: 
 tab1a = glayout(spacing=3,container=tab1) #generator
 tab1b = glayout(spacing=3,container=tab1) #create dataframe
 tab1c = glayout(spacing=20,container=tab1) #store dataframe

 refreshTab1 = function() {
  #Function for achiving entry-text: Reads into list-format:
  getEntryDataGUI <- function(asList=TRUE) {
   sname <- svalue(tab1c[1,2])
   Data <- list()
   Data[[sname]] = list(adata=list(),hdata = list())
   ln = numeric()
   for(i in 1:(length(get("popFreq",envir=mmTK))+1)) { 
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
   setEntryDataGUI(Data,sname)  
  } #end

  #load/save profiles
  f_openprof = function(h,...) {
   proffile = gfile(text="Open profile",type="open")
    if(!is.na(proffile)) {
    Data = tableReader(proffile) #load profile
    print(Data)
    Data = sample_tableToList(Data) #convert from table to list
    setEntryDataGUI(Data,names(Data)[1]) #set imported profile to GUI
   }
  }
  f_saveprof = function(h,...) {
   proffile = gfile(text="Save project",type="save")
   if(!is.na(proffile)) {
    Data = getEntryDataGUI(asList=FALSE) #get profile data as table
    Data = sample_listToTable(Data) #convert from table to list
    print(Data)
    write.table(Data,file=proffile,quote=FALSE,sep=",",row.names=FALSE) #load environment
    print(paste("Profile saved in ",proffile,sep=""))
   }
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
 f_loadkd = function(h,...) {
  loadKitList(freqpath=get("freqfolder",envir=mmTK))
  print(get("kits",envir=mmTK))
  #update combolist
 }

 #b) load/save profiles
 f_importprof = function(h,...) {
  type=h$action #get type of profile
  proffile = gfile(text=paste("Open ",type,"-profile",sep=""),type="open")
  if(!is.na(proffile)) {
   Data = tableReader(proffile) #load profile
   print(Data) #print input data
   Data = sample_tableToList(Data) #convert from table to list
   #get already stored data:
   if(type=="mix") Data2 <- getData("mix") #get data from mmTK-environment
   if(type=="ref") Data2 <- getData("ref") #get data from mmTK-environment
   if(is.null(Data2)) {
     Data2 <- Data
   } else {
     for(k in 1:length(Data)) { #for each profile
      Data2[[names(Data)[k]]] <- Data[[k]]
     }
   }
   if(type=="mix")  assign("mixData",Data2,envir=mmTK) #assign data to mmTK-environment
   if(type=="ref")  assign("refData",Data2,envir=mmTK) #assign data to mmTK-environment
  }
  #update datatable:
  if(type=="mix")  tab2b[2,1][] <- names(Data2)
  if(type=="ref")  tab2b[2,2][] <- names(Data2)
 }

 #c) 
 f_viewdata = function(h,...) {
   mixSel <- svalue(tab2b[2,1])  #get selected mixtures
   refSel <- svalue(tab2b[2,2])  #get selected references
   kit <- svalue(tab2a[2,1]) #get kit name. Must be same name as in generateEPG
   mixD <- getData(h$action) #get mixture data
   if(h$action=="mix") { #call on epg-function for each mixtures
    for(msel in mixSel) {
     subD <- mixD[[msel]]
     for(loc in names(subD$adata)) {
      if(length(grep("AMEL",loc))>0) next
      subD$adata[[loc]] <- as.numeric(subD$adata[[loc]])
     }
     plotEPG(subD,kit,msel) #plot epg's
    } 
   }
   if(h$action=="ref") { #call on epg-function for each mixtures
    refD <- getData(h$action) #get also ref-data
    for(rsel in refSel) {
     print(refD[[msel]])
    }
    #want to make a table which compares the mix-allele with ref-alleles
   }
 } 

 #type layout: 
 tab2a = glayout(spacing=5,container=tab2) #kit and population selecter
 tab2b = glayout(spacing=5,container=tab2) #evidence,ref dataframe
 tab2c = glayout(spacing=30,container=tab2) #view evidence,ref data
 tab2d = glayout(spacing=30,container=tab2) #Tasks button

 #Choose box and import button
 tab2a[1,1] = gbutton(text="Select kit directory",container=tab2a,handler=
  function(h,...) {
   dirfile = gfile(text="Select folder",type="selectdir")
   if(!is.na(dirfile)) assign("freqfolder",dirfile,envir=mmTK) #assign freqfolder
  }
 )
 tab2a[1,2] = gbutton(text="Import kit information",container=tab2a,handler=
  function(h,...) {
   loadKitList(freqpath=get("freqfolder",envir=mmTK))
   kitList <- get("kits",envir=mmTK)
   print(kitList)
   tab2a[2,1][] <- names(kitList)
  }
 )
 tab2a[2,1] <- gcombobox(items="",  selected = 0, editable = FALSE, container = tab2a, handler=
    function(h,...) {
     kitList <- get("kits",envir=mmTK)
     tab2a[2,2]['text'] <- names(kitList[[svalue(tab2a[2,1])]])
    })

 tab2a[2,2] <- gcombobox(items="",  selected = 0, editable = FALSE , container = tab2a, handler=
    function(h,...) {
     kitList <- get("kits",envir=mmTK)
     popList <- kitList[[svalue(tab2a[2,1])]][[svalue(tab2a[2,2])]] #get selected frequencies
     print(popList)
     assign("popFreq",popList,envir=mmTK) #assign popFreq get("popFreq",envir=mmTK)
    })

 #Choose box and import button
 tab2b[1,1] = gbutton(text="Import evidence",container=tab2b,handler=f_importprof,action="mix")
 tab2b[2,1] = gcheckboxgroup(items="", container = tab2b,use.table=TRUE)

 #Choose box and import button
 tab2b[1,2] = gbutton(text="Import references",container=tab2b,handler=f_importprof,action="ref")
 tab2b[2,2] = gcheckboxgroup(items="", container = tab2b,use.table=TRUE)

 #view data:
 tab2c[1,1] = gbutton(text="View evidence",container=tab2c,handler=f_viewdata,action="mix")
 tab2c[1,2] = gbutton(text="View references",container=tab2c,handler=f_viewdata,action="ref")

 #Choices:
 tab2d[1,1] = gbutton(text="Create profiles",container=tab2d,handler=
  function(h,...) {
   refreshTab1() #refresh table for creating profiles
   svalue(nb) <- 1 #change tab of notebook
  }
 )
 tab2d[1,2] = gbutton(text="Deconvolution",container=tab2d,handler=
  function(h,...) { 
   mixSel <- svalue(tab2b[2,1])  #get selected mixtures
   refSel <- svalue(tab2b[2,2])  #get selected references
   print(mixSel)
   print(refSel)
   if(length(mixSel)==0) {
    tkmessageBox(message="Please import and select mix-profile!")
   } else {
    refreshTab3(mixSel,refSel) #refresh table with selected data
    svalue(nb) <- 3 #change tab of notebook
   }
  }
 ) 

############################################################
###############Tab 3: Deconvolution setup:##################
############################################################

 #layout: 
# tab3a = glayout(spacing=1,container=tab3) #Profile selecter
 tab3b = glayout(spacing=1,container=tab3) #locus Selecter
 tab3c = glayout(spacing=1,container=tab3) #User-input selecter

 refreshTab3 = function(mixSel,refSel) {
  #mixSel,refSel is name of profiles in mixData,refData
  nM = length(mixSel) #number of mix-profiles
  nR = length(refSel) #number of ref-profiles
  mixD = getData("mix")
  refD = getData("ref") 
  locnames <- NULL
  locstart <- 3 #row of locus start - 1
  tab3b[1,1] <- glabel(text="Profile:",container=tab3b)
  tab3b[2,1] <- glabel(text="Condition:",container=tab3b)
  tab3b[3,1] <- glabel(text="",container=tab3b)
  for(nm in 1:nM) {
    tab3b[1,1 + nm] <- glabel(text=mixSel[nm],container=tab3b)
    tab3b[2,1 + nm] <- gcheckbox(text="",container=tab3b,checked=TRUE)
    subD <- mixD[[mixSel[nm]]] #select profile
    if(nm==1) locnames <- toupper(names(subD$adata)) #get loc-names
    for(i in 1:length(subD$adata)) {
     locind <- grep(names(subD$adata)[i],locnames) #get loc-index of stain:
     if(nm==1) tab3b[locstart+i,1] <- locnames[i] #insert loc-name
     tab3b[locstart+locind,1 + nm]  <- gcheckbox(text="",container=tab3b,checked=TRUE)
     #note: some data may miss some loci
    }
  }
  if(nR>0) {
   for(nr in 1:nR) {
    tab3b[1,1 + nM + nr] <- glabel(text=refSel[nr],container=tab3b)
    tab3b[2,1 + nM + nr] <- gcheckbox(text="",container=tab3b,checked=FALSE)
    subD <- refD[[refSel[nm]]] #select profile
    for(i in 1:length(subD$adata)) {
     locind <- grep(names(subD$adata)[i],locnames) #get loc-index of stain:
     if(nm==1) tab3b[locstart+i,1] <- locnames[i] #insert loc-name
     tab3b[locstart+locind,1 + nM + nr]  <- gcheckbox(text="",container=tab3b,checked=TRUE)
     #note: some data may miss some loci 
    }
   }
  } #end if refs.

  #User input info:
  tab3c[1,1] <- glabel(text="",container=tab3c)
  tab3c[2,1] <- glabel(text="Methods:",container=tab3c)
  tab3c[2,2] <- glabel(text="Models:",container=tab3c)
  tab3c[3,1] <-  gradio(items=c("Simple Search","Greedy Search","Peeloff Search"),container=tab3c,selected=3)
  tab3c[3,2] <-  gradio(items=c("Ordinary","Weighted","Generalized"),container=tab3c,selected=3)
  tab3c[3,3] <- glabel(text="        ",container=tab3c)
  tab3c[4,1] <- glabel(text="",container=tab3c)
  #options:
  tab3c[5,1] <- glabel(text="#contributors:",container=tab3c)
  tab3c[6,1] <- gedit(text="2",container=tab3c,width=3)
  tab3c[5,2] <- glabel(text="Meth. param:",container=tab3c)
  tab3c[6,2] <- gedit(text="100",container=tab3c,width=3) #startvalue
  tab3c[5,3] <- glabel(text="Threshold:",container=tab3c)
  tab3c[6,3] <- gedit(text="50",container=tab3c,width=3)
  tab3c[5,4] <- glabel(text="Allow zeroMx:",container=tab3c)
  tab3c[6,4] <- gcheckbox(text="",container=tab3c,checked=FALSE)
 
  
  tab3c[3,4] = gbutton(text="Do deconvolution!",container=tab3c,handler=
   function(h,...) { 
    #take out selected data and send them to deconvolution
    mixD2 <- list() #will take data with locnames-order
    refD2 <- list() #will take data with locnames-order
    locsel_Mix <- numeric()
    locsel_Ref <- numeric()
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
     tkmessageBox(message="Please select one mix-profile!")
     return
    }
    if(nR>0) { #note that some references may have other order of loci: This is taken into account.
     for(nr in 1:nR) { #for each selected references
      if(svalue(tab3b[2,1 + nM + nr])) { #if checked
       tmpA = list() #adata-list
       locsel = rep(FALSE,length(locnames))
       for(i in 1:length(locnames)) { #for each needed locus
         locind <- grep(locnames[i],toupper(names(refD[[refSel[nr]]]$adata))) #get loc-index of prof
         if(length(locind)==0) tmpA[[locnames[i]]] <- numeric() #assign empty if not found in mix.
         if(length(locind)==1) { #if locus of reference was found in mixture
          tmpA[[locnames[i]]] <- refD[[refSel[nr]]]$adata[[locind]] #temp-store in list
          #only accept locus if valid and checked
          if(length(tmpA[[locnames[i]]])==2 && svalue(tab3b[locstart+i,1 + nr])) locsel[i] = TRUE
         }
       } #end for each locus
       locsel_Ref = cbind(locsel_Ref , locsel) #add column
       refD2[[refSel[nr]]] = tmpA #insert list directly
      } #end if profile checked
     } #end for each ref
    } #end if any refs

    #Call for Deconvolution: Rename ref-samples first:
     rData <- NULL
     lsRef <- NULL
     condO <- NULL 
     if(length(refD2)>0) { #if any conditioned references
      rData<-refD2 
      lsRef<-locsel_Ref
      condO <- 1:length(refD2)  #conditional order of the checked references
     }
     deconvlist <- deconvolve(mixData=mixD2[[1]],nC=as.numeric(svalue(tab3c[6,1])),method=substr(svalue(tab3c[3,1]),1,4),model=substr(svalue(tab3c[3,2]),1,4),eps=as.numeric(svalue(tab3c[6,2])),locsel_Mix=locsel_Mix[,1], refData=rData,locsel_Ref=lsRef,condOrder=condO,zeroMx=svalue(tab3c[6,4]),threshT=as.numeric(svalue(tab3c[6,3])))
#     print(deconvlist)
     assign("deconvlist",deconvlist,mmTK) #store results from deconv
     #send deconvolved results to next frame
     refreshTab4(layout="Layout2") #update result-tab
     svalue(nb) <- 4 #change notetab
   } #end handle function
  ) #end button
 } #end refreshTab3

##############################################################
###############Tab 4: Deconvolution results:##################
##############################################################
 tab4a = glayout(spacing=1,container=tab4,expand=TRUE) #table layout
 tab4b = glayout(spacing=1,container=tab4,expand=TRUE) #table of results
 tab4c = glayout(spacing=1,container=tab4,expand=TRUE) #storing result

 f_savetable = function(h,...) {
   deconvlist<-get("deconvlist",envir=mmTK) #load results from environment
   if(is.null(deconvlist)) {
    tkmessageBox(message="There is no deconvolution results available.")
    return
   }
   tabfile = gfile(text="Save table",type="save")
   if(!is.na(tabfile)) {
    if(h$action=="Layout1") tab <- deconvlist$result1
    if(h$action=="Layout2") tab <- deconvlist$result2
    write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=TRUE) #load environment
    print(paste("Result table saved in ",tabfile,sep=""))
   }
  }

 refreshTab4 = function(layout="Layout2") {
   deconvlist<-get("deconvlist",envir=mmTK) #load results from environment
   if(!is.null(deconvlist)) {
    if(layout=="Layout1") restab <- deconvlist$result1
    if(layout=="Layout2") restab <- deconvlist$result2
    #gtable(restab,chosencol=1,container=TRUE)
    tab4b[1,1] <- gtable(restab,container=tab4b,multiple=TRUE)
    #tab4b[1,1] <- gtable(restab,container=tab4b,multiple=TRUE,width=500,height=500,do.autoscroll=TRUE,noRowsVisible=TRUE) #add to frame
   }
   #create table in tab4a
   tab4b[1,1] <- glabel(text="",container=tab4b)
 #  tab4b[3,2] <- gradio(items=layouts,container=tab4b,selected=grep(layout,layouts))
#   tab4b[4,2] <- gbutton(text="Update table",container=tab4b, handler=function(h,...) { refreshTab4( svalue(tab4b[3,2])  )}) #refresh table with selected layout
   tab4c[1,1] <- glabel(text="",container=tab4c)
   tab4c[2,1] <- glabel(text="             ",container=tab4c)
   tab4c[2,2] <- gbutton(text="Save table",container=tab4c,handler=f_savetable,action=svalue(tab4a[1,2]))  
 }

 layouts <- c("Layout1","Layout2")
 tab4a[1,1] <- glabel(text="Table layout:",container=tab4a)
 tab4a[1,2] <- gradio(items=layouts,container=tab4a,horisontal=TRUE, handler=function(h,...) { refreshTab4( svalue(tab4a[1,2])  )})


##############################################################
###############Tab 5: LR-calculation setup:###################
##############################################################
 


 visible(mainwin) <- TRUE


} #end funcions
