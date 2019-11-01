---
title: "Read Coffalyser Output"
author: "A Sherborne"
editor: "A Sherborne"
date: "31st October 2019"
output: html_document
---
############################################################################################################
##### This script takes in multiple digitalCofallyser output files and combines them into one          #####
##### It cleans the lab number IDs, removes references and pastes in the read count for the BRAF probe #####
############################################################################################################

#install.packages('BiocManager', repos='https://cran.ma.imperial.ac.uk/')
#BiocManager::install("xlsx", site_repository = 'https://cran.ma.imperial.ac.uk/')
#### For my Mac it required the dependencies mgcv and rJava:
#BiocManager::install("mgcv", site_repository = 'https://cran.ma.imperial.ac.uk/')
#BiocManager::install("rJava", site_repository = 'https://cran.ma.imperial.ac.uk/')

#### xlsx package seems to have some issues - one being that it jams up java. If you get an error such as:
#### Error in .jnew("org/apache/poi/xssf/usermodel/XSSFWorkbook") : 
####      java.lang.OutOfMemoryError: GC overhead limit exceeded
# run the following to clear it:
#gc()
#J("java.lang.Runtime")$getRuntime()$gc()

.libPaths() 
getwd()

library(xlsx)

############################################################################################################
##### To read just one plate and clean it, uncomment below and paste in the filename                   #####
############################################################################################################

### Read the first tab of the excel file (Final Ratios) NB. Coffalyser outputs as xls files
#Data <- read.xlsx("NamelessExperiment01.xls", 1, header=TRUE)

#Renaming sample names - looks for pattern "BP{2digits}.{2digits}_" and replaces with nothing ("") leaving the lab number
#colnames(Data) <- sub("BP\\d{2}.\\d{2}_", "", colnames(Data))

### Read the second tab of the excel file (Read Counts) to get BRAF info
#DataBRAF <- read.xlsx("NamelessExperiment01.xls", 2, header=TRUE)
#colnames(DataBRAF) <- sub("BP\\d{2}.\\d{2}_", "", colnames(DataBRAF))
#Data[Data$Gene=="BRAF",] <- DataBRAF[DataBRAF$Gene=="BRAF",]


############################################################################################################
##### Find the Coffalyser files in the working directory NB. Coffalyser outputs as xls files           #####
############################################################################################################

print(files <- list.files(pattern="*.xls$"))

Data <- data.frame()
for (i in 1:length(files)) {
    Datai <- read.xlsx(files[i], 1, header=TRUE)
    colnames(Datai) <- sub("BP\\d{2}.\\d{2}_", "", colnames(Datai))
    DataBRAFi <- read.xlsx(files[i], 2, header=TRUE)
    colnames(DataBRAFi) <- sub("BP\\d{2}.\\d{2}_", "", colnames(DataBRAFi))
    Datai[Datai$Gene=="BRAF",] <- DataBRAFi[DataBRAFi$Gene=="BRAF",]
    if (i==1) {
        Data <- Datai
    } else {
        Data <- cbind(Data, Datai[,7:ncol(Datai)])
    }
}

### Remove reference probes
DataTest <- Data[(Data$ProbeFunctions=="Test probe" | Data$ProbeFunctions=="Mutation detection probe"),]


############################################################################################################
##### Write to file. One tab with just the test probes, one with all data                              #####
############################################################################################################

currentDate <- Sys.Date()
currentDate <- format(currentDate, format="%y%m%d")
FileName <- paste(currentDate,"_Combined_dMLPA_data.xlsx",sep="")

write.xlsx2(DataTest, file=FileName, sheetName="TestProbes", col.names=TRUE, row.names=TRUE, append=FALSE)
write.xlsx2(Data, "Combined_dMLPA_data.xlsx", sheetName="AllProbes", col.names=TRUE, row.names=TRUE, append=TRUE)
