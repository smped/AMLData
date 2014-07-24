## Inspecting the AML data from DA, Mahmoud & Kyaw

## Author: Steve Pederson
## Email: stephen.pederson@adelaide.edu.au
## Commenced: 16-Apr-2014

source("loadPackages.R")

## Load & format the data
datFile <- file.path("data","data.txt")
data <- read.delim(datFile, skip=38, header=FALSE)
annot <- read.delim(datFile, nrows=38, header=FALSE, colClasses="character", stringsAsFactors=FALSE)
## Remove the first row from the annotation file as it is not useful
## Set the rownames to be the 2nd column as this defined the anonotation type
annot <- annot[-1,]
rownames(annot) <- annot[,2]
annot <- annot[,-c(1,2)]
## The annotations are set as row vectors, but R needs column vectors
Age <- as.integer(annot["Age",])
ID <- as.integer(annot["MicroarrayIDinDB",])
annot <- data.frame(ID, Age, t(annot[-c(1:2),]), stringsAsFactors=TRUE)
## Re annotate the samples so that they have an X-prefix
samples <- paste("X", substr(paste("00", ID, sep=""), start=nchar(ID), stop=nchar(ID)+2), sep="_")
annot <- data.frame(sampleID=samples, annot[,-1])
rownames(data) <- paste("ID", apply(data[,2:1], FUN=paste, MARGIN=1, collapse="_"), sep="")
data <- as.matrix(data[,-c(1:2)])
colnames(data) <- samples
n <- length(samples)
nG <- nrow(data)

## Inspect the actual data
ggData <- data.frame(sample=rep(samples, each=nrow(data)), intensity=as.vector(data))
ggplot(ggData, aes(x=intensity, colour=sample)) + 
  geom_density()+
  xlim(4, 10)
## Those expression levels seem quite low compared to most Illumina datasets

## Now for the parameters. The key ones to break the dataset down on are Source, FANC, Phenotype
FANC <- rep(1, n)
FANC[which(apply(annot[,26:28], FUN=paste, MARGIN=1, collapse="")=="NNN")] <- 0
FANC <- factor(c("wt","mut")[FANC+1], levels=c("wt","mut"))
M0 <- rep(0,n)
M0[grep("M0", annot$Phenotype)] <- 1
M1 <- rep(0,n)
M1[grep("M1", annot$Phenotype)] <- 1
M2 <- rep(0,n)
M2[grep("M2", annot$Phenotype)] <- 1
## There are no M3...
M4 <- rep(0,n)
M4[grep("M4", annot$Phenotype)] <- 1
M5 <- rep(0,n)
M5[grep("M5", annot$Phenotype)] <- 1
M_other <- rep(1, 2)
M_other <- M_other - M0 - M1 - M2 - M4 - M5
M_other[which(M_other<0)] <- 0 ## Reset any to zero that had multiple M groupings
Source <- unclass(annot$Source) - 1

## Form a model matrix
modMat <- cbind(1, unclass(FANC)-1, M0, M1, M2, M4, M5, Source)
colnames(modMat) <- c("(Intercept)", "FANCmut", "M0", "M1", "M2","M4","M5", "SourceAllG")

## Estimate the weights & fit the model
w <- arrayWeights(data, design=modMat)
png("weights.png", height=600, width=1000)
barplot(w, names.arg=samples, las=2)
abline(h=1, col="blue")
dev.off()
## 32 (X_048) & 66 (X_031) get the lowest weights
annot[c(32,66),]
fit <- lmFit(data, design=modMat, weights=w)
fit <- eBayes(fit)
