library('aroma.affymetrix')

verbose <- Arguments$getVerbose(-8, timestamp=TRUE)
name <- "MyelomaIX"
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty")
cdfs <- lapply(chipTypes, FUN=function(chipType) {AffymetrixCdfFile$byChipType(chipType)})
print(cdfs)

cesNList <- list()
chipType <- chipTypes[2]
cs <- AffymetrixCelSet$byName(name, chipType=chipType)
cs <- cs[!isDuplicated(cs)]
print(cs)

acc <- AllelicCrosstalkCalibration(cs, model="CRMAv2")
print(acc)

csC <- process(acc, verbose=verbose)
print(csC)

bpn <- BasePositionNormalization(csC, target="zero")
print(bpn)

csN <- process(bpn, verbose=verbose)
print(csN)

plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE)
print(plm)

if (length(findUnitsTodo(plm)) > 0) {
  # Fit CN probes quickly (~5-10s/array + some overhead)
  units <- fitCnProbes(plm, verbose=verbose)
  str(units)
  # int [1:945826] 935590 935591 935592 935593 935594 935595 ...

  # Fit remaining units, i.e. SNPs (~5-10min/array)
  units <- fit(plm, verbose=verbose)
  str(units)
}

ces <- getChipEffectSet(plm)
print(ces)

fln <- FragmentLengthNormalization(ces, target="zero")
print(fln)

cesNList[[chipType]] <- process(fln, verbose=verbose)
print(cesNList)


