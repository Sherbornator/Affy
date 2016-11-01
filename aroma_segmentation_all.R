### Aroma for segmentation with tumor samples normalised with normal ###

#Setting up an already preprocessed data set
rm(list=ls())
library('aroma.affymetrix')
library("aroma.cn")
library("PSCBS")
log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

name <- "MyelomaIX"
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty")
tags <- "ACC,-XY,BPN,-XY,AVG"
dsNsp <- AromaUnitTotalCnBinarySet$byName(name, tags=tags, chipType=chipTypes[1])
dsSty <- AromaUnitTotalCnBinarySet$byName(name, tags=tags, chipType=chipTypes[2])

#Set full name translator
fnt <- function(names, ...) {
  pattern <- "^([A-Z0-9]{1}[0-9]{4})[_]*(.*)"
  gsub(pattern, "\\1,\\2,\\3", names)
}
setFullNamesTranslator(dsNsp, fnt)
print(dsNsp)

setFullNamesTranslator(dsSty, fnt)
print(dsSty)



#It is for at the segmentation step you need to care about merging chip types.
#The segmentation model classes of the Aroma framework (e.g. CbsModel), will take care of the merging by simply interweaving the loci/total CN estimates from multiple chip types (if such are available for the sample currently being segmented).
#Using do[AS]CRMAv2(), you will basically get an AromaUnitTotalCnBinarySet for each chip type. 
#If you place those in an R list, e.g.

dsList <- list()
dsList[["Mapping250K_Nsp"]] <- dsNsp
dsList[["Mapping250K_Sty"]] <- dsSty
print(dsList)


#You can simply do
sm <- CbsModel(dsList)
print(sm)

fit(sm, verbose=-10)

pathname <- writeRegions(sm, verbose=verbose)

