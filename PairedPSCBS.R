### Aroma for segmentation with combined T-N pairs ###
#http://www.aroma-project.org/vignettes/PairedPSCBS-lowlevel/																			

#Setting up an already preprocessed data set
rm(list=ls())
library('aroma.affymetrix')
library("aroma.cn")
library("PSCBS")

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

name <- "MyelomaIX"
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty")
cdfs <- lapply(chipTypes, FUN=function(chipType) {AffymetrixCdfFile$byChipType(chipType)})

chipType <- chipTypes[1]

T <- c("90003", "90019", "90031", "90041", "90058", "90081", "90098", "90108", "90111", "90114", "90127", "90130", "90137", "90142", "90145", "90152", "90154", "90165", "90176", "90177", "90179", "90186", "90192", "90196", "90199", "90203", "90213", "90214", "90221", "90222", "90225", "90226", "90227", "90230", "90239", "90241", "90246", "90273", "90276", "90296", "90300", "90307", "90323", "90325", "90340", "90356", "90357", "90366", "90372", "90374", "90375", "90403", "90405", "90438", "90439", "90449", "90462", "90468", "90471", "90481", "90482", "90494", "90502", "90505", "90565", "90595", "90607", "90620", "90624", "90656", "90658", "90676", "90685", "90715", "90730", "90738", "90802", "90827", "90857", "90879", "90937", "90983", "90992", "90993", "91112", "91175", "91182", "91220", "91348", "91416", "91498", "91527", "91577", "91603", "91695", "91698", "91827", "91905", "91924")
N <- c("C0003", "C0019", "C0031", "C0041", "C0058", "C0081", "C0098", "C0108", "C0111", "C0114", "C0127", "C0130", "C0137", "C0142", "C0145", "C0152", "C0154", "C0165", "C0176", "C0177", "C0179", "C0186", "C0192", "C0196", "C0199", "C0203", "C0213", "C0214", "C0221", "C0222", "C0225", "C0226", "C0227", "C0230", "C0239", "C0241", "C0246", "C0273", "C0276", "C0296", "C0300", "C0307", "C0323", "C0325", "C0340", "C0356", "C0357", "C0366", "C0372", "C0374", "C0375", "C0403", "C0405", "C0438", "C0439", "C0449", "C0462", "C0468", "C0471", "C0481", "C0482", "C0494", "C0502", "C0505", "C0565", "C0595", "C0607", "C0620", "C0624", "C0656", "C0658", "C0676", "C0685", "C0715", "C0730", "C0738", "C0802", "C0827", "C0857", "C0879", "C0937", "C0983", "C0992", "C0993", "C1112", "C1175", "C1182", "C1220", "C1348", "C1416", "C1498", "C1527", "C1577", "C1603", "C1695", "C1698", "C1827", "C1905", "C1924")
pair = data.frame(T, N, stringsAsFactors=F)

for (i in 1:nrow(pair)){
	res <- ""
	data <- ""
	dimnames <- ""
	CT <- ""
	betaT <- ""
	betaN <- ""
	df <- ""
	fit <- ""
	segs <- ""
	chromosome <- ""
	x <- ""
	csR <- AffymetrixCelSet$byName(name, chipType=chipType)
	csR <- csR[!isDuplicated(csR)]
	csR <- csR[indexOf(csR, pair[i,])]
	print(csR)

	#apply whole allele specific CRMA in one go
	#takes a 10 minutes per pair
	res <- doASCRMAv2(csR, verbose=verbose)

	#output is to totalAndFracBData/MyelomaIX,ACC,-XY,BPN,-XY,AVG,FLN,-XY/Mapping250K_Nsp/

	data <- extractPSCNArray(res$total)
	dimnames(data)[[3]] <- names(pair)
	str(data)

	# Total CNs for the tumor relative to the matched normal
	CT <- 2 * (data[,"total","T"] / data[,"total","N"])
	
	# Allele B fractions for the tumor
	betaT <- data[,"fracB","T"]

	# Allele B fractions for the normal
	betaN <- data[,"fracB","N"]

	# Get (chromosome, position) annotation data
	ugp <- getAromaUgpFile(res$total)
	chromosome <- ugp[,1,drop=TRUE]
	x <- ugp[,2,drop=TRUE]

	# Setup data structure for Paired PSCBS
	df <- data.frame(chromosome=chromosome, x=x, CT=CT, betaT=betaT, betaN=betaN)

	#fit segments
	fit <- segmentByPairedPSCBS(df, verbose=verbose)

	segs <- getSegments(fit)
	#output that can be written to file
	print(segs)

	pathname <- write.table(segs, file = paste(pair[i,1], "Nsp_segments_pairwise_with_norm.txt", sep = ""), sep = "\t", na = "NA")
}

rm(list=ls())
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

name <- "MyelomaIX"
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty")
cdfs <- lapply(chipTypes, FUN=function(chipType) {AffymetrixCdfFile$byChipType(chipType)})

chipType <- chipTypes[2]

T <- ""
N <- ""
pair <- ""
T <- c("90003", "90019", "90031", "90041", "90058", "90081", "90098", "90108", "90111", "90114", "90127", "90130", "90136", "90137", "90142", "90145", "90152", "90154", "90165", "90176", "90177", "90179", "90186", "90196", "90199", "90203", "90212", "90213", "90214", "90221", "90222", "90225", "90226", "90227", "90230", "90239", "90241", "90245", "90246", "90273", "90276", "90280", "90296", "90300", "90307", "90323", "90325", "90340", "90356", "90357", "90366", "90372", "90374", "90375", "90403", "90405", "90438", "90439", "90449", "90462", "90468", "90471", "90481", "90482", "90494", "90502", "90505", "90565", "90595", "90607", "90624", "90656", "90658", "90676", "90685", "90715", "90730", "90738", "90827", "90857", "90879", "90937", "90983", "90992", "90993", "91112", "91175", "91182", "91190", "91220", "91348", "91416", "91498", "91527", "91577", "91603", "91695", "91698", "91827", "91890", "91905", "91924")
N <- c("C0003", "C0019", "C0031", "C0041", "C0058", "C0081", "C0098", "C0108", "C0111", "C0114", "C0127", "C0130", "C0136", "C0137", "C0142", "C0145", "C0152", "C0154", "C0165", "C0176", "C0177", "C0179", "C0186", "C0196", "C0199", "C0203", "C0212", "C0213", "C0214", "C0221", "C0222", "C0225", "C0226", "C0227", "C0230", "C0239", "C0241", "C0245", "C0246", "C0273", "C0276", "C0280", "C0296", "C0300", "C0307", "C0323", "C0325", "C0340", "C0356", "C0357", "C0366", "C0372", "C0374", "C0375", "C0403", "C0405", "C0438", "C0439", "C0449", "C0462", "C0468", "C0471", "C0481", "C0482", "C0494", "C0502", "C0505", "C0565", "C0595", "C0607", "C0624", "C0656", "C0658", "C0676", "C0685", "C0715", "C0730", "C0738", "C0827", "C0857", "C0879", "C0937", "C0983", "C0992", "C0993", "C1112", "C1175", "C1182", "C1190", "C1220", "C1348", "C1416", "C1498", "C1527", "C1577", "C1603", "C1695", "C1698", "C1827", "C1890", "C1905", "C1924")
pair = data.frame(T, N)

for (i in 1:nrow(pair)){
	res <- ""
	data <- ""
	dimnames <- ""
	CT <- ""
	betaT <- ""
	betaN <- ""
	df <- ""
	fit <- ""
	segs <- ""
	chromosome <- ""
	x <- ""
	csR <- AffymetrixCelSet$byName(name, chipType=chipType)
	csR <- csR[!isDuplicated(csR)]
	csR <- csR[indexOf(csR, pair[i,])]
	print(csR)

	#apply whole allele specific CRMA in one go
	#takes a 10 minutes per pair
	res <- doASCRMAv2(csR, verbose=verbose)

	#output is to totalAndFracBData/MyelomaIX,ACC,-XY,BPN,-XY,AVG,FLN,-XY/Mapping250K_Nsp/

	data <- extractPSCNArray(res$total)
	dimnames(data)[[3]] <- names(pair)
	str(data)

	# Total CNs for the tumor relative to the matched normal
	CT <- 2 * (data[,"total","T"] / data[,"total","N"])
	
	# Allele B fractions for the tumor
	betaT <- data[,"fracB","T"]

	# Allele B fractions for the normal
	betaN <- data[,"fracB","N"]

	# Get (chromosome, position) annotation data
	ugp <- getAromaUgpFile(res$total)
	chromosome <- ugp[,1,drop=TRUE]
	x <- ugp[,2,drop=TRUE]

	# Setup data structure for Paired PSCBS
	df <- data.frame(chromosome=chromosome, x=x, CT=CT, betaT=betaT, betaN=betaN)

	#fit segments
	fit <- segmentByPairedPSCBS(df, verbose=verbose)

	segs <- getSegments(fit)
	#output that can be written to file
	print(segs)

	pathname <- write.table(segs, file = paste(pair[i,1], "Sty_segments_pairwise_with_norm.txt", sep = ""), sep = "\t", na = "NA")
}

