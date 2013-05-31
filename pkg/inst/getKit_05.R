# TODO - Degradation method dependent degradation factor?
# TODO - Kit dependent quant-rfu factor?
# TODO - Automatically create info from GM panels etc?
# TODO - Check offset. Do not mix actual bp with min range in marker range (panel/bins).
#        Best so far: estimate the offset by taking the smallest ladder fragment i.e. 98.28.
#        and round this to an integer (98) and subtract the number of base pair for that repeat i.e. 4*9=36,
#        which gives an offset of 98-36 = 62 bp for D3. 

# New in 05: Simple test kit added.
# New in 04: ESX17 and corrected  an error in NGM.
# New in 03: interlocus.balance # Test values



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
