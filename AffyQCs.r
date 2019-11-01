#!/usr/bin/R

##### This script takes in .CEL files and cleans and normalises them using the Bioconducter package #####
##### With our dataset, two batches need to be created on upload - one with MMpprofiler arrays, the other U133s #####
##### http://www.arrayanalysis.org/ ?????

# To run in the bash terminal so R on Davros looks for user libraries:
# export R_LIBS_USER=/scratch/users/asherborne01/R/x86_64-redhat-linux-gnu-library/3.2:$R_LIBS_USER

.libPaths() 

#install.packages("affy", site_repository = 'https://cran.ma.imperial.ac.uk/')
#source("https://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("RSQLite")
#biocLite("AnnotationDbi")
#biocLite("simpleaffy")
#biocLite("sva")
library(affy)
library(simpleaffy)
library(sva)
getwd()

celsU <- read.delim("U133CelPath.txt", header = FALSE, sep = "\t")
celsM <- read.delim("MMProfilerCelPath.txt", header = FALSE, sep = "\t")
celsAddOn <- read.delim("AddOnCelPath.txt", header = FALSE, sep = "\t")

### MMProfilers ###
DataM <- ReadAffy(filenames=celsM[,1], celfile.path = './', cdfname='hgu133plus2cdf')
print("read DataM")

print(DataM)
#AffyBatch object
#size of arrays=1164x1164 features (281 kb)
#cdf=hgu133plus2cdf (54675 affyids)
#number of samples=929
#number of genes=54675
#annotation=hgu133plus2cdf
#notes=

#Checking data loaded
sampleNames(DataM)

#Renaming - looks for pattern "_@{32digits}.CEL
sampleNames(DataM) <- sub(".CEL$", "", sampleNames(DataM))

#### Removing arrays with scale factors outside 3 fold change from QC ####
badMM <- read.delim("badMMScaleFactors.txt", header = FALSE, sep = "\t")
badArrayM <- match(badMM[,1], sampleNames(DataM))
DataM <- DataM[, -badArrayM]

### Add Ons ###
DataA <- ReadAffy(filenames=celsAddOn[,1], celfile.path = './', cdfname='hgu133plus2cdf')
print("read DataA")

print(DataA)

#Checking data loaded
sampleNames(DataA)

#Renaming - looks for pattern "_@{32digits}.CEL
sampleNames(DataA) <- sub(".CEL$", "", sampleNames(DataA))


qcsM <- qc(DataM)

ppi <- 600
png(file="DataM_qcPlot_cleaning.png",width=7*ppi,height=36*ppi,res=ppi)
plot(qcsM)
dev.off()
write.table(ratios(qcsM), file="qcsM_ratios_cleaned.txt", quote=F, sep="\t", col.names=NA)
write.table(sfs(qcsM), file="qcsM_sfs_cleaned.txt", quote=F, sep="\t", col.names=NA)
write.table(percent.present(qcsM), file="qcsM_PP_cleaned.txt", quote=F, sep="\t", col.names=NA)
write.table(avbg(qcsM), file="qcsM_avbg_cleaned.txt", quote=F, sep="\t", col.names=NA)

DataM <- merge(DataM, DataA)
sampleNames(DataM)

print(DataM)

### U133s ###
DataU <- ReadAffy(filenames=celsU[,1], celfile.path = './', cdfname='hgu133plus2cdf')
print("read DataU")

print(DataU)
#AffyBatch object
#size of arrays=1164x1164 features (98 kb)
#cdf=hgu133plus2cdf (54675 affyids)
#number of samples=287
#number of genes=54675
#annotation=hgu133plus2cdf
#notes=

#Checking data loaded
sampleNames(DataU)

#Renaming - looks for pattern "_@{32digits}.CEL
sampleNames(DataU) <- sub(".CEL$", "", sampleNames(DataU))

#### Removing arrays with scale factors outside 3 fold change from QC ####
badArrayU <- match(c("17_0419", "11_0809", "11_1212", "12_0628", "15_2707", "16_1025", "16_1575", "13_2219", "11_0935", "12_2610", "14_3190", "16_0485", "16_0793", "16_1402", "16_1473", "16_1510", "17_0126", "17_0254", "17_0570", "11_0073", "11_0261", "12_0117"), sampleNames(DataU))
DataU <- DataU[, -badArrayU]

# For 329 file
# badArray <- match(c("11_0582", "11_0587", "13_1224", "13_1538", "12_0193", "12_0385"), sampleNames(Data))

#qcsU <- qc(DataU)

#ppi <- 600
#png(file="DataU_qcPlot_cleaning.png",width=7*ppi,height=14*ppi,res=ppi)
#plot(qcsU)
#dev.off()
#write.table(ratios(qcsU), file="qcsU_ratios_cleaning.txt", quote=F, sep="\t", col.names=NA)
#write.table(sfs(qcsU), file="qcsU_sfs_cleaning.txt", quote=F, sep="\t", col.names=NA)
#write.table(percent.present(qcsU), file="qcsU_PP_cleaning.txt", quote=F, sep="\t", col.names=NA)
#write.table(avbg(qcsU), file="qcsU_avbg_cleaning.txt", quote=F, sep="\t", col.names=NA)



###Combining datasets
Data <- merge(DataM, DataU)
sampleNames(Data)

print(Data)
#qcs <- qc(Data)

#ppi <- 600
#png(file="Data_qcPlot_cleaning.png",width=7*ppi,height=14*ppi,res=ppi)
#plot(qcs)
#dev.off()
#write.table(ratios(qcs), file="qcs_ratios.txt", quote=F, sep="\t", col.names=NA)
#write.table(sfs(qcs), file="qcs_sfs.txt", quote=F, sep="\t", col.names=NA)
#write.table(percent.present(qcs), file="qcs_PP.txt", quote=F, sep="\t", col.names=NA)
#write.table(avbg(qcs), file="qcs_avbg.txt", quote=F, sep="\t", col.names=NA)

##########################################################################################
###########            Normalising                                             ###########
##########################################################################################
eset <- mas5(Data)
#eset.rma <- rma(Data)

### Log 2-ing the data and writing to file ###
exprSet.nologs = exprs(eset)
write.table(exprSet.nologs, file="191023_1169Cels_eset_mas5_matrix_nologs.txt", quote=F, sep="\t", col.names=NA)

exprSet = log(exprSet.nologs, 2)
write.table(exprSet, file="191023_1169Cels_eset_mas5_matrix.txt", quote=F, sep="\t", col.names=NA)

exprSet.nologs <- as.matrix(read.table("191023_1169Cels_eset_mas5_matrix_nologs.txt", header = TRUE, sep = "\t", row.names=1, as.is=TRUE))
SampleNames <- colnames(exprSet.nologs)
SampleNames <- gsub("^X", "",  SampleNames)
colnames(exprSet.nologs) <- SampleNames

##########################################################################################
###########            Pairwise Comparisons                                    ###########
##########################################################################################

annot <- read.delim("Comb_batch_annot.txt", header = TRUE, sep = "\t")
rownames(annot) <- annot$X
annot <- AnnotatedDataFrame(annot)
all(rownames(annot)==colnames(exprSet.nologs))
#TRUE

AnnotSet <- ExpressionSet(assayData=exprSet.nologs, phenoData=annot, annotation='hgu133plus2cdf')

get.array.subset(AnnotSet,"ArrayType","MMProfiler")
results <- pairwise.comparison(AnnotSet,"ArrayType",c("MMProfiler","U133"))
print(results)

write.table(means(results), file="AnnotSet_PairComp_means.txt", quote=F, sep="\t", col.names=NA)
write.table(fc(results), file="AnnotSet_PairComp_fc.txt", quote=F, sep="\t", col.names=NA)
