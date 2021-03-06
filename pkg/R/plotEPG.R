#' @title plotEPG
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description EPG plotter maintained by Oskar Hansson.
#' @export
#' @details Plots peak height with corresponding allele for a sample.
#' @param Data List of adata- and hdata-elements.
#' @param kit name of kit.
#' @param sname label text.


#wrapper function for plotting mixture data as epg profile
#Using Oskars functions. Put generateEPG here
plotEPG <- function(Data,kitname,sname="") {
 #Data is list with allele and height data. Only one sample!
 #for selected sample:


getKit<-function(kitNameOrIndex=NULL, showMessages=TRUE) {

	# TODO: Error in getKit(typingKit, showMessages = TRUE) : 
 	# unused argument(s) (showMessages = TRUE)
	# No error with the line below:
	if(showMessages){
		print(paste("showMessages:", showMessages))
	}
	# This function returns the following information for a supported kit:
	# Short kit name.
	# Full kit name.
	# Marker/locus names.
	# Dye for each marker/locus
	# Start offset, in base pairs, for each marker.
	# Size of repeating unit, in base pairs, for each marker.

	# If no matching kit or kit index is found NA is returned.
	# If NULL a vector of available kits is printed and NA returned.

	# Available kits. Must match else if construct.
	kits<-c("SGMPlus","NGM","ESX17", "TestKit")
	
	# Check if NULL
	if (is.null(kitNameOrIndex)) {

		# Print available kits
		if (showMessages){
			print("Available kits:")
			print(kits)
		}
		kit<-NA

	# String provided.
	} else {

		# Check if number or string.
		if (is.numeric(kitNameOrIndex)) {

			# Set index to number.
			index<-kitNameOrIndex

		} else {

			# Find matching kit index (case insensitive)
			index<-match(toupper(kitNameOrIndex),toupper(kits))

		}

		# No matching kit.
		if (is.na(index)) {
			
			# Print available kits
			if (showMessages){
				print("No matching kit!")
				print("Available kits:")
				print(kits)
			}
			kit<-NA

		# Assign matching kit information.
		} else if (index == 1) {
		
			# SGM Plus
			kit<-list(
				shortName = kits[index],
				fullName = "AmpFlSTR SGM Plus PCR Amplification Kit",
				locus = c("D3S1358","vWA","D16S539","D2S1338","Amelogenin","D8S1179","D21S11","D18S51","D19S433","TH01","FGA"),
				dye = c("B","B","B","B","G","G","G","G","Y","Y","Y"),
				offset = c(65,112,213,233,100,95,90,236,70,148,146),
				repeatUnit = c(4,4,4,4,4,4,4,4,4,4,4),
				locusBalanceMean = c(1,0.7,1.1,0.9,0.8,1,0.8,0.7,1.2,0.9,0.9),
				locusBalanceSd = c(0.1,0.07,0.11,0.09,0.08,0.1,0.08,0.07,0.12,0.09,0.09)
			)
		} else if (index == 2) {

			kit<-list(
				shortName = kits[index],
				fullName = "AmpFlSTR NGM PCR Amplification Kit",
				locus = c("D10S1248","vWA","D16S539","D2S1338","Amelogenin","D8S1179","D21S11","D18S51","D22S1045","D19S433","TH01","FGA","D2S441","D3S1358","D1S1656","D12S391"),
				dye = c("B","B","B","B","G","G","G","G","Y","Y","Y","Y","R","R","R","R"),
				offset = c(72,149,223.6,281.6,100,117.9,178.8,259.5,76,122.3,176.4,221.6,74.5,114.4,170,225),
				repeatUnit = c(4,4,4,4,9,4,4,4,3,4,4,4,4,4,4),
				locusBalanceMean = c(1,0.7,1.1,0.9,0.8,1,0.8,0.7,1.2,0.9,0.9,1,1,0.6,0.8),
				locusBalanceSd = c(0.1,0.07,0.11,0.09,0.08,0.1,0.08,0.07,0.12,0.09,0.09,0.1,0.1,0.06,0.08)
			)

		} else if (index == 3) {

			kit<-list(
				shortName = kits[index],
				fullName = "PowerPlex ESX 17 System",
				locus = c("AMEL","D3S1358","TH01","D21S11","D18S51","D10S1248","D1S1656","D2S1338","D16S539","D22S1045","vWA","D8S1179","FGA","D2S441","D12S391","D19S433","SE33"),
				dye = c("B","B","B","B","B","G","G","G","G","Y","Y","Y","Y","R","R","R","R"),
				offset = c(87,103,152,203,286,83,137,197,273,79,124,203,264,88,130,193,267),
				repeatUnit = c(6,4,4,4,4,4,4,4,4,3,4,4,4,4,4,4,4),
				locusBalanceMean = c(1,0.7,1.1,0.9,0.8,1,0.8,0.7,1.2,0.9,0.9,1,1,0.6,0.8,1),
				locusBalanceSd = c(0.1,0.07,0.11,0.09,0.08,0.1,0.08,0.07,0.12,0.09,0.09,0.1,0.1,0.06,0.08,0.5)
			)

		} else if (index == 4) {
		  
		  kit<-list(
		    shortName = kits[index],
		    fullName = "A Simple Test Kit",
		    locus = c("A","B","C","D"),
		    dye = c("B","B","G","R"),
		    offset = c(72,149,223.6,281.6),
		    repeatUnit = c(2,3,4,5),
		    locusBalanceMean = c(1,0.7,1.1,0.9),
		    locusBalanceSd = c(0.1,0.07,0.11,0.09)
		  )
		  
		  # No matching index. Available kits vector does not match else if construction.
		} else {

			if (showMessages){
				print("ERROR: Kit details not defined in function!")
			}
			kit<-NA
		}	

	}

	# Return kit information (or NA)
	return(kit)

}

generateEPG<-function(alleleList, peakHeightList, dyeVector=NULL, locusVector=NULL, offsetVector=NULL,
	repeatUnitVector=NULL, sampleName=NULL, typingKit=NULL, drawBoxPlots=TRUE, drawPeaks=TRUE){

	# This function generates an electropherogram (EPG) like plot from lists of allele names and peak heights.
	# Default values are used if no typing kit is specified.
	# If a typing kit is specified some or all of its properties can be overriden.

	# INPUT PARAMETERS:

	# alleleList - a list of alleles names.
	# peakHeightList - a list of peak heights.
	# dyeVector - a vector specifying the color of markers.
	# locusVector - a vector specifying the marker names.
	# offsetVector - a vector specifying the marker start offsets in base pairs.
	# repeatUnitVector - a vector specifying the repeat unit size in base pairs.
	# sampleName - title on the EPG.
	# typingKit - kit short name or index number as specified in the function 'getKit'.

	# NB! If 'locusVector', 'dyeVector', 'offsetVector' or 'repeatUnitVector' is provided
	# they override the information in a matching 'typingKit' IF they are of the same lenght.


	# FLAGS:

	kitFound <- FALSE


	# CONSTANTS:

	# Default sample name / epg header:
	defSampleName <- "EPG"

	# Default repeat unit in base pairs:
	defaultRepeatUnit <- 4

	# Default locus name (will be suffixed with a number):
	defaultLocusName <- "Locus "

	# Default marker spacing in base pairs:
	defaultMarkerSpacing <- 100

	# Available letter abreviations for colors:
	colorLetters <- c("X", "R", "B", "G", "Y")

	# Numeric values corresponding to color abreviations:
	# NB! There are 8 colors that can be represented by a single number character, palette():
	# 1="black", 2="red", 3="green3", 4="blue", 5="cyan", 6="magenta", 7="yellow", 8="gray" 
	colorNumbers <- c(1, 2, 4, 3, 7)


	# GRAPH CONSTANTS:
	
	# Graph title:
	# sampleName is used as title.
	
	# X axis label:
	xlabel <- "base pair"
	
	# Y axis label:
	ylabel <- "peak height (rfu)"	

	# Peak half width in base pair (width of peak = 2 * peakHalfWidth):
	peakHalfWidth <- 0.6

	# Relative Allele name text size:	
	alleleNameTxtSize <- 0.7

	# Distance between the highest peak in a plot and the plot border (1.04 = 4% margin).
	yMarginTop <- 1.04


	# VERIFICATIONS
	# REQUIRED PARAMETERS:	

	# Check if allele name list is provided.
	if(is.null(alleleList) || !is.list(alleleList))	{

		stop("'alleleList' must be a list giving the alleles per marker/locus")
	}

	# Check if peak height list is provided.
	if(is.null(peakHeightList) || !is.list(peakHeightList)) {

		stop("'peakHeightList' must be a list giving the peak heights per marker/locus")
	}

	# VERIFICATIONS
	# OPTIONAL PARAMETERS:	

	# Check if sample name is provided.
	if (is.null(sampleName)) {
		sampleName <- defSampleName
	}

	# Check if kit name is provided.
	if (!is.null(typingKit)){
		
		# Get kit information.
		kit <- getKit(typingKit, showMessages=TRUE)
	
		# Check if kit was found.	
		if (is.na(kit[[1]])) {
	
			print("Kit not found.")
			print("Using default values.")
			kitFound <- FALSE
			
		} else { # Kit found. Save information in vectors.
		
			kitFound <- TRUE

			# Marker/locus names.
			locusVectorKit <- kit$locus

			# Dye for each marker/locus
			dyeVectorKit <- kit$dye

			#Base pair start offset for each marker.
			offsetVectorKit <- kit$offset

			#Size in base pair of repeating unit for each marker.
			repeatUnitVectorKit <- kit$repeatUnit
		}
	}

	# If kit is specified...
	if (kitFound){
		# ... check length of allele list.
		numberOfAlleles <- length(alleleList)
		numberOfAllelesKit <- length(locusVectorKit)
		if (numberOfAlleles != numberOfAllelesKit) {
			numberOfAlleles <- numberOfAlleles + 1
			# Extend with missing alleles (add NAs)
			for (index in numberOfAlleles:numberOfAllelesKit) {
				alleleList[[index]] <- c(NA)
			}
		}
		# ... and length of peak height list.
		numberOfPeaks <- length(peakHeightList)
		numberOfAllelesKit <- length(locusVectorKit)
		if (numberOfPeaks != numberOfAllelesKit) {
			numberOfPeaks <- numberOfPeaks + 1
			# Extend with missing peaks (add NAs)
			for (index in numberOfPeaks:numberOfAllelesKit) {
				peakHeightList[[index]] <- c(NA)
			}
		}
	}

	# Check dye vector.
	if (is.null(dyeVector)) {
		if (kitFound) {
			dyeVector <- dyeVectorKit
		} else {
			# No dye vector is specified. Use default color.
			dyeVector <- rep("X", length(alleleList))
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(dyeVector) != length(dyeVectorKit)) {
				print("'dyeVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				dyeVector <- dyeVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}

	# Check locus vector.
	if (is.null(locusVector)) {
		if (kitFound) {
			locusVector <- locusVectorKit
		} else {
			# No locus vector is specified. Use default name + number suffix.
			locusVector <- paste(rep(defaultLocusName, length(alleleList)), seq(1:length(alleleList)))
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(locusVector) != length(locusVectorKit)) {
				print("'locusVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				locusVector <- locusVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}
	# Check offset vector.
	if (is.null(offsetVector)) {
		if (kitFound) {
			offsetVector <- offsetVectorKit
		} else {
			# No offset vector is specified. Use default marker spacing.
			offsetVector <- seq(0, length(alleleList) * defaultMarkerSpacing, by = defaultMarkerSpacing)
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(offsetVector) != length(offsetVectorKit)) {
				print("'offsetVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				offsetVector <- offsetVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}
	# Check repeatUnit vector.
	if (is.null(repeatUnitVector)) {
		if (kitFound) {
			repeatUnitVector <- repeatUnitVectorKit
		} else {
			# No repeatUnit vector is specified. Use default repeat unit length.
			repeatUnitVector <- rep(defaultRepeatUnit, length(alleleList))
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(repeatUnitVector) != length(repeatUnitVector)) {
				print("'repeatUnitVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				repeatUnitVector <- repeatUnitVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}


	# PREPARE PARAMETERS

	# Convert dye vector to numeric vector.
	colorVector <- colorNumbers[match(dyeVector, colorLetters)]
	
	# Get number of colors.
	noColors <- nlevels(factor(dyeVector))

	# Get colors
	colors <- levels(factor(dyeVector))
	colors <- colorNumbers[match(colors, colorLetters)]


	# INITIATE VARIABLES

	# Create lists:
	basePairTmpLst <- list()
	basePairList <- list()
	phListByColor <- list()
	allelesByColorList <- list()


	# SORT DATA ACCORDING TO COLOR CHANNEL and
	# CONVERT ALLELE NAMES TO FRAGMENT LENGTH IN BASE PAIRS

	# Loop over all color channels.
	for (color in 1:noColors){

		# Boolean vector indicating selected markers (same color).
 		selectedMarkers <- colorVector == colors[color]

		# Extract all alleles in the same color channel.
		allelesByColor <- alleleList[selectedMarkers]

		# Extract all peak heights in the same color channel.
		peakHeightsByColor <- peakHeightList[selectedMarkers]

		# Extract all marker offsets in the same color channel.
		offsetByColor <- offsetVector[selectedMarkers]

		# Extract all repeat unit sizes in the same color channel.
		repeatUnitByColor <- repeatUnitVector[selectedMarkers]

		# Loop over all markers in that color channel.
		for (marker in 1:length(allelesByColor)){

			# Get alleles for the current marker.
			alleleValue <- allelesByColor[[marker]]
			# Check presence of X/Y.
			indexOfXY <- grep("[X,x,Y,y]", alleleValue)
			if (length(indexOfXY)) {
				alleleValue <- toupper(alleleValue)
				# Use 1 and 2 for X and Y.
				alleleValue <- sub("X", 1, alleleValue)
				alleleValue <- sub("Y", 2, alleleValue)
				alleleValue <- as.numeric(alleleValue)
			} else { #added �B
			 alleleValue <- as.numeric(alleleValue) #added �B
                  } #added �B

			# Convert all allele names in current marker to base pairs.
			basePairTmpLst[[marker]] <- offsetByColor[marker] + floor(alleleValue) * repeatUnitByColor[marker] + (alleleValue %% 1) * 10

		}

		# Add basepair to list.
		for(row in 1:length(allelesByColor)){

			# Avoid 'subscript out of bounds' error.
			if (length(basePairList)<color){
				basePairList[[color]] <- basePairTmpLst[[row]]
				phListByColor[[color]] <- peakHeightsByColor[[row]]
				allelesByColorList[[color]] <- allelesByColor[[row]]
			} else {
				basePairList[[color]] <- c(basePairList[[color]], basePairTmpLst[[row]])
				phListByColor[[color]] <- c(phListByColor[[color]], peakHeightsByColor[[row]])
				allelesByColorList[[color]] <- c(allelesByColorList[[color]], allelesByColor[[row]])
			}
		}
	}


	# CREATE GRAPH

	# Set up the plot window according to the number of color channels.
	par(mfrow = c(noColors, 1))

	# Make x axis cross at y = 0 (i.e. remove the 4% margin at both ends of the y axis).
	par(yaxs = "i")

	# Reduce the spacing between the plots.
	# c(bottom, left, top, right) is the number of lines of margin to be specified on the four sides of the plot.
	# The default is c(5, 4, 4, 2) + 0.1
	par(mar = c(3, 4, 2, 2) + 0.1)

	# Define lower and upper bound for the x axis.
	xMin <- .Machine$integer.max 
	xMax <- 0
	for (row in 1:length(basePairList)){

#print("basePairList[[row]]")
#print(basePairList[[row]])

		xChannelMin <- sapply(na.omit(basePairList[[row]]),min)
		if (length(xChannelMin) > 0) {
			xChannelMin <- min(xChannelMin)
			xVal <- is.finite(xChannelMin)
		} else {
			xVal <- FALSE
		}
		if (xMin > xChannelMin && xVal == TRUE) {
			xMin <- xChannelMin
		}
		xChannelMax <- sapply(na.omit(basePairList[[row]]),max)
		if (length(xChannelMax) > 0) {
			xChannelMax <- max(xChannelMax)
			xVal <- is.finite(xChannelMax)
		} else {
			xVal <- FALSE
		}
		if (xMax < xChannelMax && xVal == TRUE) {
			xMax <- xChannelMax
		}
	}

	#Loop over all color channels.
	for (color in 1:noColors){

		# Get alleles and peak heights for current marker.
		bpVector <- basePairList[[color]]
		phList <- phListByColor[[color]]
#print("phList")
#print(phList)
		# Create blank plot with axes.
#		yMax <- sapply(na.omit(phList),max)
		yMax <- sapply(phList[!is.na(phList)],max)
#print("yMax1")
#print(yMax)
		noData <- FALSE
		if (length(yMax) == 0) {
#print("length yMax == 0")
			yMax <- 1 # Minimal height.
			noData <- TRUE
		}
		yMax <- max(yMax)
#		yMax <- sapply(na.omit(yMax),max)
#print("Ymax2")
#print(yMax)

#		noData <- FALSE
		if (is.infinite(yMax) || is.na(yMax)) {
#print("yMax infinite or NA")
#print("yMax")
			yMax <- 1 # Minimal height.
			noData <- TRUE
		}
		plot(c(xMin, xMax), c(0, yMax), type="n", ylim = c(0, yMax * yMarginTop), ann = FALSE)
#		plot(c(xMin, xMax), c(min(phList), max(phList)), type="n", ylim = c(0, yMax * yMarginTop), ann = FALSE)


		# Write text if no data.
		if (noData) {
			text(xMax / 1.4, yMax / 2, labels="No data", cex = 1.5)
		}

		# Create a title.
		if (color == 1) {
		title(main = sampleName, col.main = "red", font.main = 4)
		}	

		# Label the y axes.
		title(ylab = ylabel)

		# Label the x axis.
		# pos values of 1, 2, 3 and 4 indicate positions below, left, above and right of the coordinate.
		mtext(paste(xlabel), side = 1, line = 2, adj = 0, cex = 0.8)

		# Write allele names under the alleles.
		# The additional par(xpd=TRUE) makes it possible to write text outside of the plot region.
		text(bpVector, 0, labels = allelesByColorList[[color]], cex = alleleNameTxtSize, pos = 1, xpd = TRUE) 

		# Loop over all peaks.
		for (peak in 1:length(bpVector)){

			# Check if boxplot is to be drawn. Do not dra if only one value.
			if (drawBoxPlots && length(phList[[peak]]) > 1) {
			# Draw boxplots showing the distribution of peak heights.
			boxplot(phList[peak], add=TRUE, at=bpVector[peak],
				border=colors[color],pars=list(boxwex=peakHalfWidth*20, axes=FALSE))
			}


			if (drawPeaks) {
			# Create corners of peak triangle.
			xCords <- c(bpVector[peak] - peakHalfWidth, bpVector[peak], bpVector[peak] + peakHalfWidth)
			yCords <- c(0, lapply(phList[peak],mean), 0)
			# Plot peaks as filled polygons.
			polygon(xCords, yCords, col = colors[color])
			}
		}
	}
 par(mfrow = c(1, 1)) #
} #end plot function

 if(all(length(unlist(Data$hdata))==0)) { 
  gmessage(message="There is no height data in sample!",title="Error",icon="error")
 } else if(all(length(unlist(Data$adata))==0)) {
  gmessage(message="There is no allele data in sample!",title="Error",icon="error")
 } else {
   #fix order prior:
   kit <- getKit(kitname, showMessages=FALSE)
   if (!is.na(kit[[1]])) { #the kit was found.
      sname <- paste0(kitname," - ",sname)
      kitlocs <- toupper(kit$locus)
      AMELname <- kitlocs[grep("AM",kitlocs,fixed=TRUE)] #get amel-name
      if(!is.null(AMELname)) { #change amel-name in dataset if AMEL was found
       names(Data$adata)[grep("AM",toupper(names(Data$adata)),fixed=TRUE)] <- AMELname
       names(Data$hdata)[grep("AM",toupper(names(Data$hdata)),fixed=TRUE)] <- AMELname
      }
      #Insert missing loci:
      missloc <- kitlocs[!kitlocs%in%names(Data$adata)] #get missing loci
      for(loc in missloc) 	Data$adata[loc] = ""
      missloc <- kitlocs[!kitlocs%in%names(Data$hdata)] #get missing loci
      for(loc in missloc) 	Data$hdata[loc] = 0
	Data$adata <- Data$adata[kitlocs] #get correct order
	Data$hdata <- Data$hdata[kitlocs] #get correct order
   }
  generateEPG(typingKit=kitname,alleleList=Data$adata,peakHeightList=Data$hdata, locusVector=names(Data$adata),sampleName=sname, drawBoxPlots=FALSE, drawPeaks=TRUE)
 }
}

