## Inspect the Raw data as it doesn't quite look right using the version sent through
##

source("loadPackages.R")
source("extraFunctions.R")

rawPrbFile <- file.path(getwd(), "data","iltp5031_Sample_Probe_Profile - Copy.txt")
file.exists(rawPrbFile)
prbData <- read.delim(rawPrbFile, skip=7, header=TRUE)
str(prbData)
head(prbData)
## Form a a matrix of the column names
cNameMat <- matrix(colnames(prbData)[-c(1:2,1347)], ncol=8, byrow=TRUE)
colnames(cNameMat) <- c("MIN_Signal","AVG_Signal","MAX_Signal","NARRAYS","ARRAY_STDEV", "BEAD_STDEV","Avg_NBEADS", "Detection")
## Check the contents
boxplot(prbData[,cNameMat[,1]] / prbData[,cNameMat[,2]])
boxplot(prbData[,cNameMat[,3]] / prbData[,cNameMat[,2]])
## The contents of MIN/AVG/MAX are the same for each signal...

## Load data. This is bead summary data
(hdr <- read.delim(rawPrbFile, nrow=5))
## It does look like this is bg corrected but not normalised
sumPrbData <- readBeadSummaryData(rawPrbFile, skip=7, # We need to skip the first 7 rows in this file & we need to specifiy the correct columns
                               columns = list(exprs = "AVG_Signal", se.exprs = "BEAD_STDEV", Detection = "Detection", nObservations = "Avg_NBEADS")) 
## The samples are in a different order now:
ID <- as.character(phenoData(sumPrbData)$sampleID)
samples <- paste("X", substr(paste("00", ID, sep=""), start=nchar(ID), stop=nchar(ID)+2), sep="_")
phenoData(sumPrbData)$sampleID <- samples
rownames(pData(sumPrbData)) <- colnames(exprs(sumPrbData)) <- samples
## That avoids any mis-referencing of the columns

## The data definiately needs quantile normalising
normPrbData <- normaliseIllumina(sumPrbData, mathod="quantile", transform="log2")

## Estimate the weights
prbWeights <- 1/se.exprs(normPrbData)^2*normPrbData@assayData$nObservations
colnames(prbWeights) <- samples

## Check the undetectable probes:
detArr <- array(0, dim=c(dim(normPrbData)[1:2], 2), dimnames=list(featureNames(normPrbData), samples, c("p","w")))
detArr[,,1] <- Detection(normPrbData)
detArr[,,2] <- prbWeights
detPVals <- apply(detArr, FUN=Stouffer, MARGIN=1)
detPrbs <- featureNames(normPrbData)[which(detPVals >0.05)]
(nG <- length(detPrbs)) # [1] 36417
## So we lost about 23%. That's probably a good result

