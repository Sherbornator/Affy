library('aroma.affymetrix')
library('aroma.cn')
#For rowMins
library('matrixStats')
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

rootPath <- 'totalAndFracBData'
rootPath <- Arguments$getReadablePath(rootPath)

dataSet <- 'MyelomaIX,ACC,-XY,BPN,-XY,AVG'
#load dataset
ds <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType='Mapping250K_Sty', paths=rootPath)

setFullNamesTranslator(ds, function(names, ...) {
  pattern <- '^([A-Z0-9]{1}[0-9]{4})[_]*(.*)'
  gsub(pattern, '\\1,\\2,\\3', names)
})
print(ds)

#Extract the normals
#N.B. Need to do this bespoke for each batch of normal samples!
normals <- c('C0003', 'C0019', 'C0031', 'C0041', 'C0058', 'C0081', 'C0098', 'C0108', 'C0111', 'C0114', 'C0127', 'C0130', 'C0136', 'C0137', 'C0142', 'C0145', 'C0152', 'C0154', 'C0165', 'C0176', 'C0177', 'C0179', 'C0186', 'C0196', 'C0199', 'C0203', 'C0212', 'C0213', 'C0214', 'C0221', 'C0222', 'C0225', 'C0226', 'C0227', 'C0230', 'C0239', 'C0241', 'C0245', 'C0246', 'C0273', 'C0276', 'C0280', 'C0296', 'C0300', 'C0307', 'C0323', 'C0325', 'C0340', 'C0356', 'C0357', 'C0366', 'C0372', 'C0374', 'C0375', 'C0403', 'C0405', 'C0438', 'C0439', 'C0449', 'C0462', 'C0468', 'C0471', 'C0481', 'C0482', 'C0494', 'C0502', 'C0505', 'C0565', 'C0595', 'C0607', 'C0624', 'C0656', 'C0658', 'C0676', 'C0685', 'C0715', 'C0730', 'C0738', 'C0827', 'C0857', 'C0879', 'C0937', 'C0983', 'C0992', 'C0993', 'C1112', 'C1175', 'C1182', 'C1190', 'C1220', 'C1348', 'C1416', 'C1498', 'C1527', 'C1577', 'C1603', 'C1695', 'C1698', 'C1827', 'C1890', 'C1905', 'C1924')
dsN <- ds[normals]
print(dsN)

#Naive genotype calling and associated confidence scores
fullname <- paste(c(getFullName(dsN), 'NGC'), collapse=',')
chipType <- getChipType(dsN, fullname=FALSE)
outPath <- file.path("callData", fullname, chipType)

units <- NULL
if (is.null(units)) {
  df <- dsN[[1]]
  units <- seq(length=nbrOfUnits(df))
}

adjust <- 1.5

# Identify units on ChrX and ChrY
ugp <- getAromaUgpFile(dsN)
units23 <- getUnitsOnChromosome(ugp, 23)
is23 <- is.element(units, units23)
units24 <- getUnitsOnChromosome(ugp, 24)
is24 <- is.element(units, units24)

for (k in 1:102){
	dfN <- dsN[[k]]

	tags <- getTags(dfN)
	tags <- setdiff(tags, 'fracB')
	tags <- c(tags, 'genotypes')
	fullname <- paste(c(getName(dfN), tags), collapse=',')

	filename <- sprintf('%s.acf', fullname)
	gcPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE)

	csTags <- c(tags, 'confidenceScores')
	fullname <- paste(c(getName(dfN), csTags), collapse=',')
	filename <- sprintf('%s.acf', fullname)
	csPathname <- Arguments$getWritablePathname(filename, path=outPath, mustNotExist=FALSE)

	if (isFile(gcPathname) && isFile(csPathname)) {
  		next
	}

	betaN <- dfN[units,1,drop=TRUE]

	# Call gender
	gender <- callXXorXY(betaN[is23], betaN[is24], adjust=adjust, from=0, to=1)

	# Call genotypes
	naValue <- as.double(NA)
	fit <- NULL
	mu <- rep(naValue, times=length(units))
	cs <- rep(naValue, times=length(units))

	if (gender == 'XY') {
  		# All but ChrX & ChrY in male
  		isDiploid <- (!(is23 | is24))
  		use <- which(isDiploid)
  		muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, from=0, to=1,
                                        verbose=less(verbose,10))
  		fit <- attr(muT, 'modelFit')
  		mu[use] <- muT
  		use <- which(!isDiploid)
  		muT <- callNaiveGenotypes(betaN[use], cn=1, adjust=adjust, from=0, to=1,
                                         verbose=less(verbose,10))
  		mu[use] <- muT
	} else {
  		# All but ChrY in female
  		isDiploid <- (!is24)
  		use <- which(isDiploid)
  		muT <- callNaiveGenotypes(betaN[use], cn=2, adjust=adjust, from=0, to=1,
                                        verbose=less(verbose,10))
  		fit <- attr(muT, 'modelFit')
  		mu[use] <- muT
	}
	print(table(mu, exclude=NULL))

	# Translate genotype calls in fracB space to (AA,AB,BB,...)
	calls <- rep(as.character(NA), times=length(mu))
	calls[mu ==   0] <- 'AA'
	calls[mu == 1/2] <- 'AB'
	calls[mu ==   1] <- 'BB'
	print(table(calls, exclude=NULL))


	# Calculate confidence scores
	a <- fit[[1]]$fitValleys$x[1]
	b <- fit[[1]]$fitValleys$x[2]
	cs[isDiploid] <- rowMins(abs(cbind(betaN[isDiploid]-a, betaN[isDiploid]-b)))
	print(table(mu, exclude=NULL))


	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Writing genotype calls (via temporary file)
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	pathname <- gcPathname
	pathnameT <- sprintf('%s.tmp', pathname)
	nbrOfUnits <- nbrOfUnits(dfN)
	gfN <- AromaUnitGenotypeCallFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits)
	footer <- readFooter(gfN)
	footer$method <- 'NaiveGenotypeCaller'
	writeFooter(gfN, footer)

	updateGenotypes(gfN, units=units, calls=calls)

	res <- file.rename(pathnameT, pathname)
	if (!isFile(pathname)) {
  		throw('Failed to rename temporary file: ', pathnameT, ' -> ', pathname)
	}
	if (isFile(pathnameT)) {
  		throw('Failed to rename temporary file: ', pathnameT, ' -> ', pathname)
	}

	gfN <- AromaUnitGenotypeCallFile(pathname)

	print(gfN)


	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	# Writing confidence scores (via temporary file)
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	pathname <- csPathname
	pathnameT <- sprintf('%s.tmp', pathname)
	nbrOfUnits <- nbrOfUnits(dfN)
	csfN <- AromaUnitSignalBinaryFile$allocate(pathnameT, platform=getPlatform(dfN), chipType=getChipType(dfN), nbrOfRows=nbrOfUnits, types='double', size=4, signed=TRUE)
	footer <- readFooter(csfN)
	footer$method <- 'NaiveGenotypeConfidenceScoreEstimator'
	writeFooter(csfN, footer)

	csfN[units, 1] <- cs

	res <- file.rename(pathnameT, pathname)
	if (!isFile(pathname)) {
  		throw('Failed to rename temporary file: ', pathnameT, ' -> ', pathname)
	}
	if (isFile(pathnameT)) {
  		throw('Failed to rename temporary file: ', pathnameT, ' -> ', pathname)
	}

	cfN <- AromaUnitSignalBinaryFile(pathname)
	print(cfN)
}
