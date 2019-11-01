#!/usr/bin/R

##### This script takes in .CEL files and looks for batch effects and normalises them using the Bioconducter package #####
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
#biocLite("limma")
library(affy)
library(simpleaffy)
library(sva)
library(limma)
getwd()


##########################################################################################
###########            Combat                                                  ###########
##########################################################################################
#covdesc <- read.delim("Comb_batch_annot.txt", header = TRUE, sep = "\t")
#x <- read.affy("covdesc")

exprSet.nologs <- as.matrix(read.table("191023_1169Cels_eset_mas5_matrix_nologs.txt", header = TRUE, sep = "\t", row.names=1, as.is=TRUE))
SampleNames <- colnames(exprSet.nologs)
SampleNames <- gsub("^X", "",  SampleNames)
colnames(exprSet.nologs) <- SampleNames

annot <- read.delim("Comb_batch_annot.txt", header = TRUE, sep = "\t")
rownames(annot) <- annot$X

#If there are any known batches, this can be given to combat
#In fact, you need to give it a batch file, even if none are known - here making them all 1s
#batch <- c(rep(1, length=890), rep(2, length=265)
#batch <- annot$ArrayType
#modcombat = model.matrix(~1, data=annot)

#Running combat NB. we want to specify a non-parametric model here (takes longer to run)
#combat_exprSet.nologs = ComBat(dat=exprSet.nologs, batch=batch, par.prior=FALSE, prior.plots=FALSE)

#Write to file
#write.table(combat_exprSet.nologs, file="191023_1169Cels_eset_mas5_combat_matrix.txt", quote=F, sep="\t", col.names=NA)


#There are lots of other useful options at the foot of page 6 of the manual: https://www.bioconductor.org/packages/devel/bioc/vignettes/sva/inst/doc/sva.pdf

mod <- model.matrix(~as.factor(ArrayType), data=annot)
#Specify the null model
mod0 <- model.matrix(~1,data=annot)

#First it identifies the number of latent factors that need to be estimated. If the sva function is called without
#the n.sv argument specified, the number of factors will be estimated for you.
#The number of factors can also be estimated using the num.sv
#n.sv <- num.sv(exprSet.nologs,mod,method="leek")
#n.sv
# The above gave a mun.sv of 1152 which broke the model!
#> svobj <- sva(exprSet.nologs,mod,mod0,n.sv=n.sv)
#Number of significant surrogate variables is:  1165 
#Iteration (out of 5 ):1  Error in solve.default(t(mod) %*% mod) : 
#  system is computationally singular: reciprocal condition number = 5.23302e-20

#svobj <- sva(exprSet.nologs,mod,mod0,n.sv=n.sv)
#Try running without the n.sv prior

svobj <- sva(exprSet.nologs,mod,mod0)
svobj

pValues <- f.pvalue(exprSet.nologs,mod,mod0)
qValues <- p.adjust(pValues,method="BH")

modSv = cbind(mod,svobj$sv)
mod0Sv = cbind(mod0,svobj$sv)
pValuesSv = f.pvalue(exprSet.nologs,modSv,mod0Sv)
qValuesSv = p.adjust(pValuesSv,method="BH")

#Write to file
#write.table(pValuesSv, file="191029_1169Cels_eset_mas5_combat_matrix_pValuesSv.txt", quote=F, sep="\t", col.names=NA)
#write.table(qValuesSv, file="191029_1169Cels_eset_mas5_combat_matrix_qValuesSv.txt", quote=F, sep="\t", col.names=NA)

##### Adjust with limma 
fit <- lmFit(exprSet.nologs,modSv)


#pheno <- read.delim("Comb_batch_annot.txt", header = TRUE, sep = "\t")
#pheno = pData(bladderEset)
