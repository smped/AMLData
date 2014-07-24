## Repeat of the previous analysis but using the updated annotations as sent through on 30-05-2014

source("loadPackages.R")
source("extraFunctions.R")

## Check the new annotations sent through
newFile <- list.files(file.path("rawData"), full.names=TRUE)[2]
newAnnot <- read.delim(newFile)
## OK, that's hard to work with. Transpose the info in Office, save to the folder data & reopen 

newFile <- file.path("data", "NewAnnotation.csv")
file.exists(newFile)
newAnnot <- read.delim(newFile)

## Quickly compare the two
colnames(annot)
colnames(newAnnot)
## There are some fields we can remove to clean the file up a little
newAnnot <- newAnnot[,-which(colnames(newAnnot) %in% c("Cytogenetics", "Genetic.Group.Risk.Stratification..ELN."))]
## We also need to remove those with an Unknown Source as there is no FANC, C1-5, MBD, Seq.*, PraimaryAML, PrimaryMDS, PrimaryMisc, Age or Sex info
newAnnot <- newAnnot[-which(newAnnot$Source=="Unknown"),]
newAnnot$Source <- newAnnot$Source[,drop=TRUE]
## Remove the Unknown Simple Karyotypes. we lose one sample with a FANC mutation by doing this.
newAnnot <- newAnnot[-which(newAnnot$SimpKaryo=="Unknown"),]
newAnnot$SimpKaryo <- newAnnot$SimpKaryo[,drop=TRUE]
## There are 3 remaining with Sex == 0, but removing those won't change anything else. 
## They are all SimpKaryo==N, Adelaide, PrimaryAML & one s Phenotype==Misc so they can help with the baseline. 
## Sex is not incorporated into the model so it's OK to leave this with conflicting data
## Some of the fields need to have any factor levels removed for which there are no longer any samples
newAnnot$Phenotype <- newAnnot$Phenotype[, drop=TRUE]
## Trim the whitespace from that KLF5 annotation
newAnnot$KLF5..Simple.[grep("HIGH", newAnnot$KLF5..Simple.)] <- "HIGH"
newAnnot$KLF5..Simple. <- newAnnot$KLF5..Simple.[, drop=TRUE]
## Have a look:
lapply(newAnnot, FUN=table)
## Get rid of the "NN" entries from PrimaryMDS & PrimaryMisc
newAnnot$PrimaryMDS[which(newAnnot$PrimaryMDS=="NN")] <- "N"
newAnnot$PrimaryMDS <- newAnnot$PrimaryMDS[, drop=TRUE]
newAnnot$PrimaryMisc[which(newAnnot$PrimaryMisc=="NN")] <- "N"
newAnnot$PrimaryMisc <- newAnnot$PrimaryMisc[, drop=TRUE]
## Have a look again:
lapply(newAnnot, FUN=table)
## That leaves how many samples?
(n <- nrow(newAnnot)) # [1] 139

## Re-label the MicroarrayIDinDB field to be a character
newAnnot$MicroarrayIDinDB <- paste(substr(rep("X_00",n), start=1, stop=5-nchar(newAnnot$MicroarrayIDinDB)), newAnnot$MicroarrayIDinDB, sep="")
useArr <- newAnnot[,"MicroarrayIDinDB"]

####################################################
## Reform all the parameters for the model matrix ##
####################################################

## Source
Source <- rep(0, n)
Source[which(newAnnot$Source=="ALLG")] <- 1 # Sets Adelaide/Unknown as the baseline

## The FAB groupings. As there is only 1 individual in each of the M6 & M7 categories, they have been dropped
M <- matrix(0, nrow=n, ncol=6, dimnames=list(c(), c("M0", "M1", "M2", "M3", "M4", "M5")))
M[grep("M0", newAnnot$Phenotype), "M0"] <- 1
M[grep("M1", newAnnot$Phenotype), "M1"] <- 1
M[grep("M2", newAnnot$Phenotype), "M2"] <- 1
M[grep("M3", newAnnot$Phenotype), "M3"] <- 1
M[grep("M4", newAnnot$Phenotype), "M4"] <- 1
M[grep("M5", newAnnot$Phenotype), "M5"] <- 1

SimpKaryo <- matrix(0, nrow=n, ncol=2, dimnames=list(c(), c("SimpKaryo.A","SimpKaryo.C"))) # N is the baseline
SimpKaryo[which(newAnnot$SimpKaryo=="A"),1] <- 1L
SimpKaryo[which(newAnnot$SimpKaryo=="C"),2] <- 1L

## Primary AML
PrimaryAML <- rep(0, n)
PrimaryAML[which(newAnnot$PrimaryAML=="Y")] <- 1

## Primary MDS
PrimaryMDS <- rep(0, n)
PrimaryMDS[which(newAnnot$PrimaryMDS=="Y")] <- 1

## Primary Misc
PrimaryMisc <- rep(0, n)
PrimaryMisc[which(newAnnot$PrimaryMisc=="Y")] <- 1

## For FLT3, NPM1, CEBPA & many other factors from here, Unknown is treated as WT
FLT3 <- matrix(0, nrow=n, ncol=2, dimnames=list(c(), c("FLT3.ITD","FLT3.TKD")))
FLT3[which(newAnnot$FLT3.ITD=="Y"),1] <- 1
FLT3[which(newAnnot$FLT3.TKD=="Y"),2] <- 1

NPM1 <- rep(0, n)
NPM1[which(newAnnot$NPM1=="Y")] <- 1

CEBPA <- rep(0, n)
CEBPA[which(!newAnnot$CEBPA %in% c("WT", "Unknown"))] <- 1

GADD45 <- matrix(0, ncol=2, nrow=n, dimnames=list(c(), c("GADD45.Low","GADD45.High")))
GADD45[grep("low", newAnnot$GADD45..Simple.),"GADD45.Low"] <- 1
GADD45[grep("high", newAnnot$GADD45..Simple.),"GADD45.High"] <- 1

KLF5 <- matrix(0, nrow=n, ncol=2)
colnames(KLF5) <- c("KLF5.Low", "KLF5.High")
KLF5[grep("HIGH",newAnnot$KLF5..Simple.), "KLF5.High"] <- 1
KLF5[grep("LOW",newAnnot$KLF5..Simple.), "KLF5.Low"] <- 1

Seq <- matrix(0, nrow=n, ncol=7)
colnames(Seq) <- paste("Seq", c("DNMT3a", "FLT3.TKD", "IDH1", "IDH2", "NPM1", "NRAS", "WT1"), sep=".")
for (i in 1:7) Seq[which(newAnnot[,colnames(Seq)[i]]=="Mutant"), i] <- 1

FANC <- matrix(0, ncol=3, nrow=n, dimnames=list(c(),c("FANC.Full.Network", "FANC.Core.Anchor.ID", "FANC.ALL")))
FANC[which(newAnnot$FANC.Full.Network == "Y"),"FANC.Full.Network"] <- 1
FANC[which(newAnnot$FANC.Core.Anchor.ID =="Y"), "FANC.Core.Anchor.ID"] <- 1
FANC[which(newAnnot$FANC.ALL == "Y"), "FANC.ALL"] <- 1
FANC.Merged <- rowMaxs(FANC) # Group the sparate mutations together

MBD <- matrix(0, nrow=n, ncol=3)
colnames(MBD) <- paste("MBD", c(".Full.Network", ".s.and.Likes", "1.4.MECP2"), sep="")
for (i in 1:3) MBD[which(newAnnot[,colnames(MBD)[i]]=="Y"),i] <- 1

METC <- rep(0, n)
METC[which(newAnnot$METC=="Y")] <- 1

## There are only 1, 1 & 2 samples from C2, C3 & C4 respectively so these columns are a little redundant
C1 <- C5 <- rep(0, n)
C1[which(newAnnot$C1=="Y")] <- 1
C5[which(newAnnot$C5=="Y")] <- 1

fullMat <- cbind(Source, M, SimpKaryo, PrimaryAML, PrimaryMDS, PrimaryMisc, FLT3, NPM1, CEBPA, GADD45, KLF5, Seq, FANC[,1:2], MBD, METC, C1, C5)
mcor <- cor(fullMat)
png("correlations2.png", height=800, width=800)
corrplot(mcor)
dev.off()

## The inverse correlations between Source & GADD45/KLF5 juust mean that these are Unknowns for all the ALLG samples
lapply(newAnnot[which(newAnnot$Source=="ALLG"), c("GADD45..Simple.", "KLF5..Simple.")], FUN=table)
lapply(newAnnot[which(!newAnnot$Source=="ALLG"), c("GADD45..Simple.", "KLF5..Simple.")], FUN=table)
## METC is still very correlated with C1
mcor["METC","C1"]

## We now have only 12 N sample for Primary AML & 127 Ys
lapply(newAnnot[which(newAnnot$PrimaryAML=="N"), c("SimpKaryo", "Phenotype", "PrimaryMDS", "PrimaryMisc", "FANC.Full.Network", "FANC.Core.Anchor.ID", paste("C", 1:5, sep=""))], FUN=table)
## 3 of these are SimpKaryo==N, so we may get some confounding. Leave it out of the model
## All of the Ns have been reclassified to PrimaryMDS or Primary Misc. If we include these but leave Primary AML out of the model that should work

##########################
## Perform the analysis ##
##########################
modMat <- cbind(1, Source, M, SimpKaryo, PrimaryMDS, PrimaryMisc, FANC[,1], C1, C5) ## This is one of the most informative but simple models
colnames(modMat)[c(1, 13)] <- c("Intercept", "FFN")
rownames(modMat) <- as.character(newAnnot$MicroarrayIDinDB)

## Estimate the weights 
arrWeights <- arrayWeights(exprs(normPrbData)[detPrbs,useArr], design=modMat)# , weights=prbWeights[detPrbs,useArr]) # Don't use the probe-level weights!!!
barplot(arrWeights, names.arg=useArr, las=2)
abline(h=1, col="blue", lty=2)
wtMat <- arrWeights #* prbWeights[detPrbs,useArr]

## Fit the model
fit <- lmFit(exprs(normPrbData)[detPrbs,useArr], design=modMat, weights=wtMat)
fit <- eBayes(fit)
## get the top FANC genes
topTabFFN <- topTable(fit, coef="FFN", n=nG)
length(which(topTabFFN$adj.P.Val < 0.05)) # [1] 285
topTabFFN <- topTabFFN[order(topTabFFN$t, decreasing=TRUE),]
gNames <- sapply(rownames(topTabFFN), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
head(topTabFFN)
tail(topTabFFN)

## Export this file for Zeya sorted on the t-statistic
write.csv(cbind(gNames, topTabFFN), "FFN.Genelist2.csv")

## Now for the actial interaction analysis, take the top genes with an FDR = 0.1 & fit a more detailed model
## The potential interactions should be between the SimpKaryo & FFN categories
sigFFN <- rownames(topTabFFN)[which(topTabFFN$adj.P.Val<0.1)]
FFN.lm <- vector("list", length(unique(unlist(sigFFN)))) # From objects to hold the model co-efficients & the r-squared values
names(FFN.lm) <- sigFFN
for (gn in sigFFN) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0 + FFN:SimpKaryo.A + FFN:SimpKaryo.C ) # Set the initial fully parameterised model, with only the FFN interactions
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  FFN.lm[[gn]] <- step(tempLm) # Keep the final model using an automated step-wise model selection process
}

## Get the coefficients & adjust the p-values for this subset of coefficients
FFNcoefs <- array(NA, dim=c(length(sigFFN), 5, 3), dimnames=list(sigFFN, c("Estimate", "Std.Error","t", "p", "adj.P"), c("FFN", "FFN:SimpKaryo.A", "FFN:SimpKaryo.C")))
FFNcoefs[,1:4,"FFN"] <- t(sapply(FFN.lm, FUN=extractCoefs.lm, coef="FFN"))
FFNcoefs[,1:4,"FFN:SimpKaryo.A"] <- t(sapply(FFN.lm, FUN=extractCoefs.lm, coef="SimpKaryo.A:FFN"))
FFNcoefs[,1:4,"FFN:SimpKaryo.C"] <- t(sapply(FFN.lm, FUN=extractCoefs.lm, coef="SimpKaryo.C:FFN"))
FFNcoefs[,5,] <- p.adjust(FFNcoefs[,4,], method="fdr")
## Have a look to see what happens if those with a rawP >0.05 are excluded
FFNcoefs[,"Estimate",][which(FFNcoefs[,"p",]>0.05)] <- 0

## Get predicted values of logFC for each gene
## What we're looking for is the logFC due to FANC for each subtype so we'll need to have the baseline FANC (SimpKaryo.N)
## with the other two columns being this baseline + the interaction terms:
FFNpred <- cbind(FFNcoefs[,"Estimate","FFN"], 
                 FFNcoefs[,"Estimate","FFN"] + FFNcoefs[,"Estimate","FFN:SimpKaryo.A"],
                 FFNcoefs[,"Estimate","FFN"] + FFNcoefs[,"Estimate","FFN:SimpKaryo.C"])
colnames(FFNpred) <- paste("SimpKaryo", c("N","A","C"), sep=".")
rownames(FFNpred) <- sigFFN
dim(FFNpred) # [1] 700   3
## Remove the irrelevant genes
FFNpred <- FFNpred[-which(rowMins(FFNcoefs[,"adj.P",])>0.05),] # Those with no significant effects
FFNpred <- FFNpred[-which(rowMaxs(abs(FFNpred)) < log2(2)),] # Those with no large fold-change
dim(FFNpred) # [1] 61   3
## Now sort the genes based onsome sort of clustering to makea mice heatmap
FFNpred <- FFNpred[hclust(dist(FFNpred))$order,]

## Setup a heatmap
Gene <- sapply(rownames(FFNpred), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
FFNggData <- data.frame(Probe = factor(rownames(FFNpred), levels=rownames(FFNpred)),
                        Gene = factor(Gene, levels=unique(Gene)),
                        logFC = as.vector(FFNpred),
                        SimpleKaryotype = factor(rep(c("Normal","Abnormal","Complex"), each=nrow(FFNpred)), levels=c("Normal", "Abnormal","Complex")))
sz=4
ggplot(FFNggData, aes(x=SimpleKaryotype, y=Probe, fill=logFC)) +
  ggtitle("Estimated log2 fold-change for:\nFull Network Mutations Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(names(Gene), Gene, sep=" / "), expand=c(0,0)) +
  theme(axis.text.x=element_text(size=sz), axis.text.y=element_text(size=sz-1.5),
        axis.title.x=element_text(size=sz), axis.title.y=element_text(size=sz),
        legend.title=element_text(size=sz), legend.text=element_text(size=sz),
        plot.title=element_text(size=sz+1), axis.ticks=element_line(size=0.15),
        legend.position="bottom", legend.title.align=0,
        legend.key.width=unit(.7,"cm"), legend.key.height=unit(0.15,"cm"))
ggsave("FFNbyKaryotype2.png", height=15, width=8, units="cm")
