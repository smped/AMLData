## Look for interactions between FANC & any of the other annotations

source("loadPackages.R")
source("extraFunctions.R")

## Remove any arrays with Unknown for the Karyotype
annot <- annot[-which(annot$SimpKaryo=="Unknown"),]
useArr <- annot[,1]
n <- length(useArr)

## Reform the model matrix for the extra removals:
Source <- rep(0, n)
Source[which(annot$Source=="ALLG")] <- 1 # Sets Adelaide/Unknown as the baseline
## Work through each of the other annotations from here...

## The FAB groupings. As there is only 1 individual in each of the M6 & M7 categories, they have been dropped
M <- matrix(0, nrow=n, ncol=6, dimnames=list(c(), c("M0", "M1", "M2", "M3", "M4", "M5")))
M[grep("M0", annot$Phenotype), "M0"] <- 1
M[grep("M1", annot$Phenotype), "M1"] <- 1
M[grep("M2", annot$Phenotype), "M2"] <- 1
M[grep("M3", annot$Phenotype), "M3"] <- 1
M[grep("M4", annot$Phenotype), "M4"] <- 1
M[grep("M5", annot$Phenotype), "M5"] <- 1

PrimaryAML <- rep(0, n)
PrimaryAML[which(annot$PrimaryAML=="Y")] <- 1

FLT3 <- matrix(0, nrow=n, ncol=2, dimnames=list(c(), c("FLT3.ITD","FLT3.TKD")))
FLT3[which(annot$FLT3.ITD=="Y"),1] <- 1
FLT3[which(annot$FLT3.TKD=="Y"),2] <- 1
FLT3.Merged <- rowMaxs(FLT3)

NPM1 <- rep(0, n)
NPM1[which(annot$NPM1=="Y")] <- 1

GADD45 <- matrix(0, ncol=2, nrow=n, dimnames=list(c(), c("GADD45.Low","GADD45.High")))
GADD45[grep("low", annot$GADD45..Simple),"GADD45.Low"] <- 1
GADD45[grep("high", annot$GADD45..Simple),"GADD45.High"] <- 1

FANC <- matrix(0, ncol=3, nrow=n, dimnames=list(c(),c("FANC.Full.Network", "FANC.Core.Anchor.ID", "FANC.ALL")))
FANC[which(annot$FANC.Full.Network == "Y"),"FANC.Full.Network"] <- 1
FANC[which(annot$FANC.Core.Anchor.ID =="Y"), "FANC.Core.Anchor.ID"] <- 1
FANC[which(annot$FANC.ALL == "Y"), "FANC.ALL"] <- 1
FANC.Merged <- rowMaxs(FANC) # Group the sparate mutations together

MBD <- rep(0, n)
MBD[grep("Y",apply(annot[,28:30], FUN=paste, MARGIN=1, collapse=""))] <- 1

METC <- rep(0, n)
METC[grep("Y", annot$METC)] <- 1

## There are only 1, 1 & 2 samples from C2, C3 & C4 so these columns are a little redundant
C1 <- C5 <- rep(0, n)
C1[which(annot$C1=="Y")] <- 1
C5[which(annot$C5=="Y")] <- 1

SimpKaryo <- matrix(0, nrow=n, ncol=2, dimnames=list(c(), c("SimpKaryo.A","SimpKaryo.C"))) # N is the baseline
SimpKaryo[which(annot$SimpKaryo=="A"),1] <- 1L
SimpKaryo[which(annot$SimpKaryo=="C"),2] <- 1L

## Specify the a relatively simple model
modMat <- cbind(1, Source, M, SimpKaryo, PrimaryAML, FANC[,1], C1, C5) ## This is one of the most informative but simple models
colnames(modMat)[c(1, 12)] <- c("Intercept", "FFN")
rownames(modMat) <- useArr

## Estimate the weights 
arrWeights <- arrayWeights(exprs(normPrbData)[detPrbs,useArr], design=modMat, weights=prbWeights[detPrbs,useArr])
barplot(arrWeights, names.arg=useArr, las=2)
abline(h=1, col="blue", lty=2)
wtMat <- arrWeights * prbWeights[detPrbs,useArr]

## Fit the model
fit <- lmFit(exprs(normPrbData)[detPrbs,useArr], design=modMat, weights=wtMat)
fit <- eBayes(fit)
## get the top FANC genes
topTabFFN <- topTable(fit, coef="FFN", n=nG)
length(which(topTabFFN$adj.P.Val < 0.05)) # [1] 238
topTabFFN <- topTabFFN[order(topTabFFN$t, decreasing=TRUE),]
gNames <- sapply(rownames(topTabFFN), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))

## Export this file for Zeya sorted on the t-statistic
write.csv(cbind(gNames, topTabFFN), "FFN.Genelist.csv")

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

## Get predicted values of logFC for each gene
## What we're looking for is the logFC due to FANC for each subtype so we'll need to have the baseline FANC (SimpKaryo.N)
## with the other two columns being this baseline + the interaction terms:
FFNpred <- cbind(FFNcoefs[,"Estimate","FFN"], 
                 FFNcoefs[,"Estimate","FFN"] + FFNcoefs[,"Estimate","FFN:SimpKaryo.A"],
                 FFNcoefs[,"Estimate","FFN"] + FFNcoefs[,"Estimate","FFN:SimpKaryo.C"])
colnames(FFNpred) <- paste("SimpKaryo", c("N","A","C"), sep=".")
rownames(FFNpred) <- sigFFN
dim(FFNpred) # [1] 619   3
## Remove the irrelevant genes
FFNpred <- FFNpred[-which(rowMins(FFNcoefs[,"adj.P",])>0.05),] # Those with no significant effects
FFNpred <- FFNpred[-which(rowMaxs(abs(FFNpred)) < log2(1.5)),] # Those with no large fold-change
dim(FFNpred) # [1] 140   3
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
ggsave("FFNbyKaryotype.png", height=15, width=8, units="cm")

## As an alternative, just fit the full interaction model
modMat2 <- cbind(1, Source, M, SimpKaryo, PrimaryAML, FANC[,1], C1, C5, FANC[,1] * SimpKaryo)
colnames(modMat2)[c(1, 12, 15, 16)] <- c("Intercept", "FFN", "FFN:SimpKaryo.A", "FFN:SimpKaryo.C")
## Fit the model
fit2 <- lmFit(exprs(normPrbData)[detPrbs,useArr], design=modMat2, weights=wtMat)
fit2 <- eBayes(fit2)
sigFFN <- rownames(coef(fit2))[unique(c(which(topTable(fit2, coef="FFN", n=nG, sort.by="none")$adj.P.Val < 0.1),
                                                which(topTable(fit2, coef="FFN:SimpKaryo.A", n=nG, sort.by="none")$adj.P.Val < 0.1),
                                                which(topTable(fit2, coef="FFN:SimpKaryo.C", n=nG, sort.by="none")$adj.P.Val < 0.1)))]
length(sigFFN) # [1] 689

## Now fit these genes removing all non-significant terms & repeat the above
FFN.lm <- vector("list", length(unique(unlist(sigFFN)))) # From objects to hold the model co-efficients & the r-squared values
names(FFN.lm) <- sigFFN
for (gn in sigFFN) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat2[,1:14]) # Form a temporary data object with the probe data & the full model matrix
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

## Get predicted values of logFC for each gene
## What we're looking for is the logFC due to FANC for each subtype so we'll need to have the baseline FANC (SimpKaryo.N)
## with the other two columns being this baseline + the interaction terms:
FFNpred <- cbind(FFNcoefs[,"Estimate","FFN"], 
                 FFNcoefs[,"Estimate","FFN"] + FFNcoefs[,"Estimate","FFN:SimpKaryo.A"],
                 FFNcoefs[,"Estimate","FFN"] + FFNcoefs[,"Estimate","FFN:SimpKaryo.C"])
colnames(FFNpred) <- paste("SimpKaryo", c("N","A","C"), sep=".")
rownames(FFNpred) <- sigFFN
dim(FFNpred) # [1] 689   3
## Remove the irrelevant genes
FFNpred <- FFNpred[which(rowMins(FFNcoefs[,"adj.P",])<0.05),] # Only those with significant effects
FFNpred <- FFNpred[which(rowMaxs(abs(FFNpred)) > log2(1.5)),] # Only those with large fold-change
dim(FFNpred) # [1] 122   3
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

## That seems to be the best approach really & those 122 are a good starting point for knowing which 
## are affected across karyotypes
## To export the csv, we need co-efficients, adjusted p-values & logFC estimates
write.csv(cbind(Gene, FFNcoefs[names(Gene), "Estimate",], FFNcoefs[names(Gene), "adj.P", ], FFNpred), "FFNByKaryotype.csv")

## Redo the above but being more restrictive on log fold change (>1)
temp <- which(rowMaxs(abs(FFNpred))>1)
length(temp) # [1] 48
## Setup a heatmap
FFNggData <- data.frame(Probe = factor(rownames(FFNpred)[temp], levels=rownames(FFNpred)[temp]),
                        Gene = factor(Gene[temp], levels=unique(Gene[temp])),
                        logFC = as.vector(FFNpred[temp,]),
                        SimpleKaryotype = factor(rep(c("Normal","Abnormal","Complex"), each=length(temp)), levels=c("Normal", "Abnormal","Complex")))
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
ggsave("FFNbyKaryotype3.png", height=10, width=8, units="cm")



#################################################################################################################################
#################################################################################################################################
## Repeat the above but breaking the interactions down by FAB grouping

## As an alternative, just fit the full interaction model
modMatFAB <- cbind(1, Source, M, SimpKaryo, PrimaryAML, FANC[,1], C1, C5, FANC[,1] * M[,c(2,3,5,6)])
colnames(modMatFAB)[c(1, 12, 15, 16, 17, 18)] <- c("Intercept", "FFN", "FFN:M1", "FFN:M2", "FFN:M4", "FFN:M5")
## Fit the model
fitFAB <- lmFit(exprs(normPrbData)[detPrbs,useArr], design=modMatFAB, weights=wtMat)
fitFAB <- eBayes(fitFAB)
sigFFN.FAB <- rownames(coef(fitFAB))[unique(c(which(topTable(fitFAB, coef="FFN", n=nG, sort.by="none")$adj.P.Val < 0.1),
                                        which(topTable(fitFAB, coef="FFN:M1", n=nG, sort.by="none")$adj.P.Val < 0.1),
                                        which(topTable(fitFAB, coef="FFN:M2", n=nG, sort.by="none")$adj.P.Val < 0.1),
                                        which(topTable(fitFAB, coef="FFN:M4", n=nG, sort.by="none")$adj.P.Val < 0.1),
                                        which(topTable(fitFAB, coef="FFN:M5", n=nG, sort.by="none")$adj.P.Val < 0.1)))]
length(sigFFN.FAB) # [1] 1189

## Now fit these genes removing all non-significant terms & repeat the above
FFN.FAB.lm <- vector("list", length(unique(unlist(sigFFN.FAB)))) # From objects to hold the model co-efficients & the r-squared values
names(FFN.FAB.lm) <- sigFFN.FAB
for (gn in sigFFN.FAB) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMatFAB[,1:14]) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0 + FFN:M1 + FFN:M2 + FFN:M4 + FFN:M5) # Set the initial fully parameterised model, with only the FFN interactions
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  FFN.FAB.lm[[gn]] <- step(tempLm) # Keep the final model using an automated step-wise model selection process
}

## Get the coefficients & adjust the p-values for this subset of coefficients
FFN.FABcoefs <- array(NA, dim=c(length(sigFFN.FAB), 5, 5), dimnames=list(sigFFN.FAB, c("Estimate", "Std.Error","t", "p", "adj.P"), c("FFN", "FFN:M1", "FFN:M2", "FFN:M4", "FFN:M5")))
FFN.FABcoefs[,1:4,"FFN"] <- t(sapply(FFN.FAB.lm, FUN=extractCoefs.lm, coef="FFN"))
FFN.FABcoefs[,1:4,"FFN:M1"] <- t(sapply(FFN.FAB.lm, FUN=extractCoefs.lm, coef="M1:FFN"))
FFN.FABcoefs[,1:4,"FFN:M2"] <- t(sapply(FFN.FAB.lm, FUN=extractCoefs.lm, coef="M2:FFN"))
FFN.FABcoefs[,1:4,"FFN:M4"] <- t(sapply(FFN.FAB.lm, FUN=extractCoefs.lm, coef="M4:FFN"))
FFN.FABcoefs[,1:4,"FFN:M5"] <- t(sapply(FFN.FAB.lm, FUN=extractCoefs.lm, coef="M5:FFN"))
FFN.FABcoefs[,5,] <- p.adjust(FFN.FABcoefs[,4,], method="fdr")

## Get predicted values of logFC for each gene
## What we're looking for is the logFC due to FANC for each subtype so we'll need to have the baseline FANC (SimpKaryo.N)
## with the other two columns being this baseline + the interaction terms:
FFN.FABpred <- cbind(FFN.FABcoefs[,"Estimate","FFN"], 
                 FFN.FABcoefs[,"Estimate","FFN"] + FFN.FABcoefs[,"Estimate","FFN:M1"],
                 FFN.FABcoefs[,"Estimate","FFN"] + FFN.FABcoefs[,"Estimate","FFN:M2"],
                 FFN.FABcoefs[,"Estimate","FFN"] + FFN.FABcoefs[,"Estimate","FFN:M4"],
                 FFN.FABcoefs[,"Estimate","FFN"] + FFN.FABcoefs[,"Estimate","FFN:M5"])
colnames(FFN.FABpred) <- paste("M", c(1,2,4,5), sep="")
rownames(FFN.FABpred) <- sigFFN.FAB
dim(FFN.FABpred) # [1] 1189   8
## Remove the irrelevant genes
FFN.FABpred <- FFN.FABpred[which(rowMins(FFN.FABcoefs[,"adj.P",])<0.05),] # Only those with significant effects
FFN.FABpred <- FFN.FABpred[which(rowMaxs(abs(FFN.FABpred)) > 1),] # Only those with large fold-change
dim(FFN.FABpred) # [1] 77   5
## Now sort the genes based on some sort of clustering to makea mice heatmap
FFN.FABpred <- FFN.FABpred[hclust(dist(FFN.FABpred))$order,]
## Setup a heatmap
Gene <- sapply(rownames(FFN.FABpred), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
FFN.FABggData <- data.frame(Probe = factor(rownames(FFN.FABpred), levels=rownames(FFN.FABpred)),
                        Gene = factor(Gene, levels=unique(Gene)),
                        logFC = as.vector(FFN.FABpred),
                        SimpleKaryotype = factor(rep(c("Misc", "M1", "M2", "M4", "M5"), each=nrow(FFN.FABpred)), levels=c("Misc", "M1", "M2", "M4", "M5")))
sz=4
ggplot(FFN.FABggData, aes(x=SimpleKaryotype, y=Probe, fill=logFC)) +
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
ggsave("FFNbyFAB.png", height=15, width=8, units="cm")

## 
## To export the csv, we need co-efficients, adjusted p-values & logFC estimates
write.csv(cbind(Gene, FFN.FABcoefs[names(Gene), "Estimate",], FFN.FABcoefs[names(Gene), "adj.P", ], FFN.FABpred), "FFNByFAB.csv")


#################################################################################################################################
#################################################################################################################################

## Repeat for the Core.Anchor.ID mutants
## First we need to remove the FFN mutants which are not FCA mutants:
r <- setdiff(which(FANC[,1]==1), which(FANC[,2]==1))
modMatFCA <- cbind(1, Source, M, SimpKaryo, PrimaryAML, FANC[,2], C1, C5, FANC[,2]*SimpKaryo)[-r,]
colnames(modMatFCA)[c(1, 12, 15, 16)] <- c("Intercept", "FCA", "SimpKaryo.A:FCA", "SimpKaryo.C:FCA")
rownames(modMatFCA) <- useArr[-r]
colSums(modMatFCA)
## That would only leave 3 samples for the Normal Karyotype so we can't really check for interactions this time.
modMatFCA <- modMatFCA[,1:14]

fitFCA <- lmFit(exprs(normPrbData)[detPrbs,useArr[-r]], design=modMatFCA, weights=wtMat[,-r])
fitFCA <- eBayes(fitFCA)
sigFCA <- rownames(coef(fitFCA))[which(topTable(fitFCA, coef="FCA", n=nG, sort.by="none")$adj.P.Val < 0.05)]
length(sigFCA) # [1] 246

## Now remove any extraneous terms from the model
## Now fit these genes removing all non-significant terms & repeat the above
FCA.lm <- vector("list", length(unique(unlist(sigFCA)))) # From objects to hold the model co-efficients & the r-squared values
names(FCA.lm) <- sigFCA
for (gn in sigFCA) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr[-r]], modMatFCA) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0) # Set the initial fully parameterised model, with only the FFN interactions
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr[-r]]) # Evaluate the initial model
  FCA.lm[[gn]] <- step(tempLm) # Keep the final model using an automated step-wise model selection process
}

## Get the coefficients & adjust the p-values for this subset of coefficients
FCAcoefs <- matrix(NA, nrow=length(sigFCA), ncol=5, dimnames=list(sigFCA, c("Estimate", "Std.Error","t", "p", "adj.P")))
FCAcoefs[,1:4] <- t(sapply(FCA.lm, FUN=extractCoefs.lm, coef="FCA"))
FCAcoefs[,5] <- p.adjust(FCAcoefs[,4], method="fdr")
## Just keep the big movers:
FCAcoefs <- FCAcoefs[which(abs(FCAcoefs[,1])>log2(1.5)),]
FCAcoefs <- FCAcoefs[order(FCAcoefs[,1], decreasing=FALSE),]
## Do a barplot
barplot(FCAcoefs[,1], las=2, horiz=TRUE)
GeneFCA <- sapply(rownames(FCAcoefs), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
ggDataFCA <- data.frame(Probe = factor(rownames(FCAcoefs), levels=rownames(FCAcoefs)),
                        Gene = factor(GeneFCA, levels=unique(GeneFCA)),
                        logFC = FCAcoefs[,"Estimate"],
                        sd = FCAcoefs[,"Std.Error"])
sz=4
ggplot(ggDataFCA, aes(y=logFC, x=Probe, fill=logFC)) +
  geom_bar(position="dodge", stat="identity") + 
  ggtitle("Estimated log2 fold-change for \nCore/Anchor/ID Mutant Vs WT")+ 
  geom_errorbar(aes(ymin = logFC-sd, ymax=logFC+sd, width=0.5), position=position_dodge(width=0.9), colour="grey80", size=0.2) +
  coord_flip() + 
  scale_fill_gradient2(low="blue", high="red") +
  scale_x_discrete(labels=paste(names(GeneFCA), GeneFCA, sep=" / "), expand=c(0,0)) +
  theme(axis.text = element_text(colour = "grey50", size = sz*0.8),
        axis.line = element_line(colour = "grey80", size=0.5),
        axis.title.y = element_blank(),
        axis.ticks = element_line(colour = "grey80", size=0.15),
        axis.title.x = element_text(colour = "grey50", size=sz),
        title = element_text(face = "bold"),
        panel.background = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(colour="grey80", size=0.25),
        legend.title=element_text(size=sz, colour = "grey50"), legend.text=element_text(size=sz, colour="grey50"),
        plot.title=element_text(size=sz+1), 
        legend.position="bottom", legend.title.align=1,
        legend.key.width=unit(.7,"cm"), legend.key.height=unit(0.15,"cm")) 
ggsave("FCA.png", height=12, width=7, units="cm")

topTabFCA <- topTable(fitFCA, n=nG, coef="FCA")
topTabFCA <- topTabFCA[order(topTabFCA$t, decreasing=TRUE),]
gNames <- sapply(rownames(topTabFCA), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
write.csv(cbind(gNames, topTabFCA), "FCA.Genelist.csv")








###################################################################################################################################
###################################################################################################################################
##                                            Initial ideas below............                                                    ##
###################################################################################################################################
###################################################################################################################################

## Form the model matrix using 3-way interactions for FANC only
modMat <- cbind(1, Source, M, SimpKaryo, PrimaryAML, FLT3, NPM1, GADD45, FANC[,1:2], MBD, METC, C1, C5, 
                FANC[,1]*SimpKaryo, FANC[,2]*SimpKaryo, FANC[,1:2]*PrimaryAML, SimpKaryo*PrimaryAML,
                FANC[,1]*SimpKaryo*PrimaryAML, FANC[,2]*SimpKaryo*PrimaryAML)
colnames(modMat)[c(1, 23:ncol(modMat))] <- c("Intercept","FFN:SimpKaryo.A", "FFN:SimpKaryo.C", "FCA:SimpKaryo.A", "FCA:SimpKaryo.C", 
                                             "FFN:PrimaryAML", "FCA:PrimaryAML", "SimpKaryo.A:PrimaryAML", "SimpKaryo.C:PrimaryAML",
                                             "FFN:SimpKaryo.A:PrimaryAML", "FFN:SimpKaryo.C:PrimaryAML", "FCA:SimpKaryo.A:PrimaryAML", "FCA:SimpKaryo.C:PrimaryAML")
colnames(modMat)[grep("FANC", colnames(modMat))] <- c("FFN", "FCA") # Use these abbreviations for consistency with the interaction terms
rownames(modMat) <- useArr
colSums(modMat)
## That looks about as complex as it can be made whilst still retaining meaning.
## Potential probelms are that of the 7*C5 samples, 5 are PrimaryAML, leaving only 2 to establish the baseline
## Similarly for FCA, as there are only 3 samples to establish a baseline for FCA:SimpKaryo.N

## Reform the weights
## Estimate the weights 
arrWeights <- arrayWeights(exprs(normPrbData)[detPrbs,useArr], design=modMat, weights=prbWeights[detPrbs,useArr])
barplot(arrWeights, names.arg=useArr, las=2)
abline(h=1, col="blue", lty=2)
wtMat <- arrWeights * prbWeights[detPrbs,useArr]

## Fit the model
fit <- lmFit(exprs(normPrbData)[detPrbs,useArr], design=modMat, weights=wtMat)
fit <- eBayes(fit)

## Check the FANC.Full.Network probes & just collect the significant probes
## for all FFN comparisons without any further info
coefFFN <- grep("FFN", colnames(modMat))
sigFANC.FFN <- vector("list", length(coefFFN)) #
names(sigFANC.FFN) <- colnames(modMat)[coefFFN]
for (i in coefFFN) {
  coef <- colnames(modMat)[i]
  topTab <- topTable(fit, coef=coef, number=nG)
  sigFANC.FFN[[coef]] <- rownames(topTab)[which(topTab$adj.P.Val<0.05)]
}
length(unique(unlist(sigFANC.FFN))) # [1] 57


#########################################################
## Now do the step-wise model selection for each probe ##
#########################################################
modMat2 <- modMat[,c(1:22)] # Without interactions! Add them to the fm object to ensure that the step process keeps them "linked"

## FANC.Full.Network...
FANC.FFN.lm <- vector("list", length(unique(unlist(sigFANC.FFN)))) # From objects to hold the model co-efficients & the r-squared values
names(FANC.FFN.lm) <- unique(unlist(sigFANC.FFN))
for (gn in unique(unlist(sigFANC.FFN))) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat2) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0
               + FFN:SimpKaryo.A + FFN:SimpKaryo.C + FFN:PrimaryAML
               + FCA:SimpKaryo.A + FCA:SimpKaryo.C + FCA:PrimaryAML
               + PrimaryAML:SimpKaryo.A + PrimaryAML:SimpKaryo.C
               + FFN:PrimaryAML:SimpKaryo.A + FFN:PrimaryAML:SimpKaryo.C
               + FCA:PrimaryAML:SimpKaryo.A + FCA:PrimaryAML:SimpKaryo.C) # Set the initial fully parameterised model, with only the FFN interactions
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  FANC.FFN.lm[[gn]] <- step(tempLm) # Keep the final model using an automated step-wise model selection process
}

nT <- 6 # The number of FFN related terms in the models
FANC.FFN.out <- array(0, dim=c(length(FANC.FFN.lm), nT , 5), 
                      dimnames=list(names(FANC.FFN.lm),
                                    c("FFN", "SimpKaryo.A:FFN", "SimpKaryo.C:FFN", "PrimaryAML:FFN", "SimpKaryo.A:PrimaryAML:FFN", "SimpKaryo.C:PrimaryAML:FFN"),
                                    c("Estimate", "Std.Error", "t", "P.Value","adj.P.Value")))
for (i in 1:nT) FANC.FFN.out[,i,1:4] <- (t(sapply(FANC.FFN.lm, FUN=extractCoefs.lm, coef=dimnames(FANC.FFN.out)[[2]][i])))
FANC.FFN.out[,,"adj.P.Value"] <- p.adjust(FANC.FFN.out[,,"P.Value"], method="fdr")
rmv <- which(apply(FANC.FFN.out[,,"adj.P.Value"], FUN=min, MARGIN=1)>0.05) # Remove any withnon-significant p-values

## Form the newData for prediction. This will give 3 FANC.WT non-PrimaryAML & PrimaryAML samples for each Karyotype & replicates for FANC.FFN.mut
newData <- data.frame(Intercept=0, Source=0, M0=0, M1=0, M2=0, M3=0, M4=0, M5=0, 
                      FLT3.ITD = 0, FLT3.TKD = 0, NPM1 = 0, GADD45.Low = 0, GADD45.High = 0, FCA = 0, METC=0, MBD = 0, C1 = 0, C5 = 0,
                      SimpKaryo.A=c(0,1,0), SimpKaryo.C=c(0,0,1), PrimaryAML=c(0,0,0,1,1,1), FFN=c(rep(0,6), rep(1,6)))
pred.FANC.FFN <- t(sapply(FANC.FFN.lm, FUN=predict, newdata=newData))
pred.FANC.FFN <- pred.FANC.FFN[,7:12] - pred.FANC.FFN[,1:6] 
pred.FANC.FFN <- pred.FANC.FFN[which(rowMaxs(abs(pred.FANC.FFN))>log2(1.5)),] # remove any for which the max logFC < 1
pred.FANC.FFN <- pred.FANC.FFN[hclust(dist(pred.FANC.FFN))$order,]
Gene <- sapply(rownames(pred.FANC.FFN), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
FANC.FFN.ggData <- data.frame(Probe = factor(rownames(pred.FANC.FFN), levels=rownames(pred.FANC.FFN)),
                              Gene = factor(Gene, levels=unique(Gene)),
                              logFC = as.vector(pred.FANC.FFN),
                              SimpleKaryotype = factor(rep(c("Normal","Abnormal","Complex"), each=nrow(pred.FANC.FFN)), levels=c("Normal", "Abnormal","Complex")),
                              PrimaryAML = c(rep("non-PrimaryAML", times=nrow(pred.FANC.FFN)*3), rep("PrimaryAML", times=nrow(pred.FANC.FFN)*3)))
ggplot(FANC.FFN.ggData, aes(x=SimpleKaryotype, y=Probe, fill=logFC)) +
  facet_wrap(~PrimaryAML) +
  ggtitle("Estimated log2 fold-change for FANC.Full.Network Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(names(Gene), Gene, sep=" / "), expand=c(0,0)) 


## None of these genes are the same as the ones that turned up without the 3-way interactions!!!
## That doesn't feel right




FANC.FFN.logFC.est <- cbind(FANC.FFN.out[,"FFN","Estimate"], # logFC for Normal Karyotype
                            FANC.FFN.out[,"FFN","Estimate"] + FANC.FFN.out[,"SimpKaryo.A:FFN","Estimate"], # logFC for Abnormal.FANCmut va AbnormalWT
                            FANC.FFN.out[,"FFN","Estimate"] + FANC.FFN.out[,"SimpKaryo.C:FFN","Estimate"], # logFC for Complex.FANCmut vs ComplexWT
                            FANC.FFN.out[,"FFN","Estimate"] + FANC.FFN.out[,"PrimaryAML:FFN","Estimate"])
colnames(FANC.FFN.logFC.est) <- c("Normal", "Abnormal", "Complex")
FANC.FFN.logFC.est <- FANC.FFN.logFC.est[which(rowMaxs(abs(FANC.FFN.logFC.est))>0.5),] # Get rid of any where the maximum FC < sqrt(2) in all comparisons
FANC.FFN.logFC.est <- FANC.FFN.logFC.est[hclust(dist(FANC.FFN.logFC.est))$order,] # Cluster using hclust
Gene <- sapply(rownames(FANC.FFN.logFC.est), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
FANC.FFN.ggData <- data.frame(Probe = factor(rownames(FANC.FFN.logFC.est), levels=rownames(FANC.FFN.logFC.est)),
                              Gene = factor(Gene, levels=unique(Gene)),
                              logFC = as.vector(FANC.FFN.logFC.est),
                              SimpleKaryotype = factor(rep(colnames(FANC.FFN.logFC.est), each=nrow(FANC.FFN.logFC.est)), levels=c("Normal", "Abnormal","Complex")))
FANC.FFN.ggData$Probe <- factor(FANC.FFN.ggData$Probe, levels=FANC.FFN.ggData$Probe) # Ensure the order is not changed by ggplot
FANC.FFN.ggData$Gene <- factor(FANC.FFN.ggData$Gene, levels=FANC.FFN.ggData$Gene) # Ensure the order is not changed by ggplot
png("FFNbyKaryotype.png", height=900, width=900)
  ggplot(FANC.FFN.ggData, aes(x=SimpleKaryotype, y=Probe, fill=logFC)) +
  ggtitle("Estimated log2 fold-change for FANC.Full.Network Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(rownames(FANC.FFN.logFC.est), Gene, sep=" / ")) # + 
  theme(axis.text.y=element_text(size=3), axis.text.x=element_text(size=4))
dev.off()
ggsave("FFNbyKaryotype.png", height=10, width=10, units="cm")




#####################################################################################################################
#####################################################################################################################



## FANC.Core.Anchor.ID
FANC.FCA.lm <- vector("list", length(unique(unlist(sigFANC.FCA)))) # From objects to hold the model co-efficients & the r-squared values
names(FANC.FCA.lm) <- unique(unlist(sigFANC.FCA))
for (gn in unique(unlist(sigFANC.FCA))) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat2) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~.+0
               +FANC.Full.Network:SimpKaryo.A
               +FANC.Full.Network:SimpKaryo.C
               +FANC.Core.Anchor.ID:SimpKaryo.A
               +FANC.Core.Anchor.ID:SimpKaryo.C) # Set the initial fully parameterised model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  FANC.FCA.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
FANC.FCA.out <- array(0, 
                      dim=c(length(FANC.FCA.lm),5 , 5), 
                      dimnames=list(names(FANC.FCA.lm),
                                    c("SimpKaryo.A", "SimpKaryo.C","FANC.Core.Anchor.ID", "SimpKaryo.A:FANC.Core.Anchor.ID", "SimpKaryo.C:FANC.Core.Anchor.ID"),
                                    c("Estimate", "Std.Error", "t", "P.Value","adj.P.Value")))
FANC.FCA.out[, 1, 1:4] <- t(sapply(FANC.FCA.lm, FUN=extractCoefs, coef="SimpKaryo.A"))
FANC.FCA.out[, 2, 1:4] <- t(sapply(FANC.FCA.lm, FUN=extractCoefs, coef="SimpKaryo.C"))
FANC.FCA.out[, 3, 1:4] <- t(sapply(FANC.FCA.lm, FUN=extractCoefs, coef="FANC.Core.Anchor.ID"))
FANC.FCA.out[, 4, 1:4] <- t(sapply(FANC.FCA.lm, FUN=extractCoefs, coef="SimpKaryo.A:FANC.Core.Anchor.ID"))
FANC.FCA.out[, 5, 1:4] <- t(sapply(FANC.FCA.lm, FUN=extractCoefs, coef="SimpKaryo.C:FANC.Core.Anchor.ID"))
FANC.FCA.out[,,"adj.P.Value"] <- p.adjust(FANC.FCA.out[,,"P.Value"], method="fdr")
for (i in 1:5) {
  FANC.FCA.out[which(FANC.FCA.out[,i,"adj.P.Value"]>0.05),i,"Estimate"] <- 0
}
FANC.FCA.logFC.est <- cbind(FANC.FCA.out[,"FANC.Core.Anchor.ID","Estimate"], # logFC for Normal Karyotype
                            FANC.FCA.out[,"FANC.Core.Anchor.ID","Estimate"] + FANC.FCA.out[,4,"Estimate"], # logFC for Abnormal.FANCmut va AbnormalWT
                            FANC.FCA.out[,"FANC.Core.Anchor.ID","Estimate"] + FANC.FCA.out[,5,"Estimate"]) # logFC for Complex.FANCmut vs ComplexWT
colnames(FANC.FCA.logFC.est) <- c("Normal", "Abnormal", "Complex")
FANC.FCA.logFC.est <- FANC.FCA.logFC.est[which(rowMaxs(abs(FANC.FCA.logFC.est))>0.5),] # Get rid of any where the maximum FC < sqrt(2) in all comparisons
FANC.FCA.logFC.est <- FANC.FCA.logFC.est[hclust(dist(FANC.FCA.logFC.est))$order,] # Cluster using hclust
Gene <- sapply(rownames(FANC.FCA.logFC.est), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
FANC.FCA.ggData <- data.frame(Probe = factor(rownames(FANC.FCA.logFC.est), levels=rownames(FANC.FCA.logFC.est)),
                              Gene = factor(Gene, levels=unique(Gene)),
                              logFC = as.vector(FANC.FCA.logFC.est),
                              SimpleKaryotype = factor(rep(colnames(FANC.FCA.logFC.est), each=nrow(FANC.FCA.logFC.est)), levels=c("Normal", "Abnormal","Complex")))
FANC.FCA.ggData$Probe <- factor(FANC.FCA.ggData$Probe, levels=FANC.FCA.ggData$Probe) # Ensure the order is not changed by ggplot
FANC.FCA.ggData$Gene <- factor(FANC.FCA.ggData$Gene, levels=FANC.FCA.ggData$Gene) # Ensure the order is not changed by ggplot
png("FCAbyKaryotype.png", height=900, width=900)
ggplot(FANC.FCA.ggData, aes(x=SimpleKaryotype, y=Probe, fill=logFC)) +
  ggtitle("Estimated log2 fold-change for FANC.Core.Anchor.ID Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(rownames(FANC.FCA.logFC.est), Gene, sep=" / ")) 
dev.off()

## C1
C1.lm <- vector("list", length(unique(unlist(sigC1)))) # From objects to hold the model co-efficients & the r-squared values
names(C1.lm) <- unique(unlist(sigC1))
for (gn in unique(unlist(sigC1))) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat2) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0 + C1:M1 + C1:M2 + C1:M4 + C1:M5 + C5:M1 + C5:M2 + C5:M4) # Set the initial fully parameterised model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  C1.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
C1.out <- array(0, 
                      dim=c(length(C1.lm), 9, 5), 
                      dimnames=list(names(C1.lm),
                                    c("M1", "M2","M4", "M5", "C1", "M1:C1", "M2:C1", "M4:C1", "M5:C1"),
                                    c("Estimate", "Std.Error", "t", "P.Value","adj.P.Value")))
for (i in 1:9) {C1.out[, i, 1:4] <- t(sapply(C1.lm, FUN=extractCoefs, coef=dimnames(C1.out)[[2]][i]))}
C1.out[,,"adj.P.Value"] <- p.adjust(C1.out[,,"P.Value"], method="fdr")
for (i in 1:9) {C1.out[which(C1.out[,i,"adj.P.Value"]>0.05),i,"Estimate"] <- 0}

C1.logFC.est <- cbind(C1.out[,"C1","Estimate"], # logFC for the heterogeneous population (M0, M3, M6, M7, Unknown)
                      C1.out[,"C1","Estimate"] + C1.out[,"M1:C1","Estimate"], # logFC for M1 C1 Vs not C1
                      C1.out[,"C1","Estimate"] + C1.out[,"M2:C1","Estimate"], # logFC for M2 C1 Vs not C1
                      C1.out[,"C1","Estimate"] + C1.out[,"M4:C1","Estimate"],
                      C1.out[,"C1","Estimate"] + C1.out[,"M5:C1","Estimate"])  
colnames(C1.logFC.est) <- c("Het", paste("M", c(1,2,4,5), sep=""))
C1.logFC.est <- C1.logFC.est[which(rowMaxs(abs(C1.logFC.est))>0.5),] # Get rid of any where the maximum FC < sqrt(2) in all comparisons
C1.logFC.est <- C1.logFC.est[hclust(dist(C1.logFC.est))$order,] # Cluster using hclust
Gene <- sapply(rownames(C1.logFC.est), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
C1.ggData <- data.frame(Probe = factor(rownames(C1.logFC.est), levels=rownames(C1.logFC.est)),
                              Gene = factor(Gene, levels=unique(Gene)),
                              logFC = as.vector(C1.logFC.est),
                              FAB_Grouping = factor(rep(colnames(C1.logFC.est), each=nrow(C1.logFC.est)), levels=colnames(C1.logFC.est)))
C1.ggData$Probe <- factor(C1.ggData$Probe, levels=C1.ggData$Probe) # Ensure the order is not changed by ggplot
C1.ggData$Gene <- factor(C1.ggData$Gene, levels=C1.ggData$Gene) # Ensure the order is not changed by ggplot
png("C1byFAB.png", height=600, width=900)
ggplot(C1.ggData, aes(x=FAB_Grouping, y=Probe, fill=logFC)) +
  ggtitle("Estimated log2 fold-change for C1 Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(rownames(C1.logFC.est), Gene, sep=" / ")) 
dev.off()

## C5
C5.lm <- vector("list", length(unique(unlist(sigC5)))) # From objects to hold the model co-efficients & the r-squared values
names(C5.lm) <- unique(unlist(sigC5))
for (gn in unique(unlist(sigC5))) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat2) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0 + C5:M1 + C5:M2 + C5:M4 + C5:M5 + C5:M1 + C5:M2 + C5:M4) # Set the initial fully parameterised model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  C5.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
C5.out <- array(0, 
                dim=c(length(C5.lm), 7, 5), 
                dimnames=list(names(C5.lm),
                              c("M1", "M2","M4", "C5", "M1:C5", "M2:C5", "M4:C5"),
                              c("Estimate", "Std.Error", "t", "P.Value","adj.P.Value")))
for (i in 1:7) {C5.out[, i, 1:4] <- t(sapply(C5.lm, FUN=extractCoefs, coef=dimnames(C5.out)[[2]][i]))}
C5.out[,,"adj.P.Value"] <- p.adjust(C5.out[,,"P.Value"], method="fdr")
for (i in 1:7) {C5.out[which(C5.out[,i,"adj.P.Value"]>0.05),i,"Estimate"] <- 0}

C5.logFC.est <- cbind(C5.out[,"C5","Estimate"], # logFC for the heterogeneous population (M0, M3, M6, M7, Unknown)
                      C5.out[,"C5","Estimate"] + C5.out[,"M1:C5","Estimate"], # logFC for M1 C5 Vs not C5
                      C5.out[,"C5","Estimate"] + C5.out[,"M2:C5","Estimate"], # logFC for M2 C5 Vs not C5
                      C5.out[,"C5","Estimate"] + C5.out[,"M4:C5","Estimate"])  
colnames(C5.logFC.est) <- c("Het", paste("M", c(1,2,4), sep=""))
C5.logFC.est <- C5.logFC.est[which(rowMaxs(abs(C5.logFC.est))>0.5),] # Get rid of any where the maximum FC < sqrt(2) in all comparisons
C5.logFC.est <- C5.logFC.est[hclust(dist(C5.logFC.est))$order,] # Cluster using hclust
Gene <- sapply(rownames(C5.logFC.est), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
C5.ggData <- data.frame(Probe = factor(rownames(C5.logFC.est), levels=rownames(C5.logFC.est)),
                        Gene = factor(Gene, levels=unique(Gene)),
                        logFC = as.vector(C5.logFC.est),
                        FAB_Grouping = factor(rep(colnames(C5.logFC.est), each=nrow(C5.logFC.est)), levels=colnames(C5.logFC.est)))
C5.ggData$Probe <- factor(C5.ggData$Probe, levels=C5.ggData$Probe) # Ensure the order is not changed by ggplot
C5.ggData$Gene <- factor(C5.ggData$Gene, levels=C5.ggData$Gene) # Ensure the order is not changed by ggplot
png("C5byFAB.png", height=600, width=900)
ggplot(C5.ggData, aes(x=FAB_Grouping, y=Probe, fill=logFC)) +
  ggtitle("Estimated log2 fold-change for C5 Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(rownames(C5.logFC.est), Gene, sep=" / ")) 
dev.off()

## There are only two samples of C5 for each FAB group so this really isn't that great.
## Have a look by Karyotype...
C1.lm <- vector("list", length(unique(unlist(sigC1)))) # From objects to hold the model co-efficients & the r-squared values
names(C1.lm) <- unique(unlist(sigC1))
for (gn in unique(unlist(sigC1))) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat2) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0 + C1:SimpKaryo.A + C1:SimpKaryo.C) # Set the initial fully parameterised model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  C1.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
C1.out <- array(0, 
                dim=c(length(C1.lm), 3, 5), 
                dimnames=list(names(C1.lm),
                              c("C1", "SimpKaryo.A:C1", "SimpKaryo.C:C1"),
                              c("Estimate", "Std.Error", "t", "P.Value","adj.P.Value")))
for (i in 1:3) {C1.out[, i, 1:4] <- t(sapply(C1.lm, FUN=extractCoefs, coef=dimnames(C1.out)[[2]][i]))}
C1.out[,,"adj.P.Value"] <- p.adjust(C1.out[,,"P.Value"], method="fdr")
for (i in 1:3) {C1.out[which(C1.out[,i,"adj.P.Value"]>0.05),i,"Estimate"] <- 0}

C1.logFC.est <- cbind(C1.out[,"C1","Estimate"], # logFC for the heterogeneous population (M0, M3, M6, M7, Unknown)
                      C1.out[,"C1","Estimate"] + C1.out[,"SimpKaryo.A:C1","Estimate"], # logFC for M1 C1 Vs not C1
                      C1.out[,"C1","Estimate"] + C1.out[,"SimpKaryo.C:C1","Estimate"]) # logFC for M2 C1 Vs not C1
colnames(C1.logFC.est) <- c("Normal", "Abnormal", "Complex")
C1.logFC.est <- C1.logFC.est[which(rowMaxs(abs(C1.logFC.est))>0.5),] # Get rid of any where the maximum FC < sqrt(2) in all comparisons
C1.logFC.est <- C1.logFC.est[hclust(dist(C1.logFC.est))$order,] # Cluster using hclust
Gene <- sapply(rownames(C1.logFC.est), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
C1.ggData <- data.frame(Probe = factor(rownames(C1.logFC.est), levels=rownames(C1.logFC.est)),
                        Gene = factor(Gene, levels=unique(Gene)),
                        logFC = as.vector(C1.logFC.est),
                        SimpleKaryotype = factor(rep(colnames(C1.logFC.est), each=nrow(C1.logFC.est)), levels=c("Normal", "Abnormal","Complex")))
C1.ggData$Probe <- factor(C1.ggData$Probe, levels=C1.ggData$Probe) # Ensure the order is not changed by ggplot
C1.ggData$Gene <- factor(C1.ggData$Gene, levels=C1.ggData$Gene) # Ensure the order is not changed by ggplot
png("C1byFAB.png", height=600, width=900)
ggplot(C1.ggData, aes(x=SimpleKaryotype, y=Probe, fill=logFC)) +
  ggtitle("Estimated log2 fold-change for C1 Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(rownames(C1.logFC.est), Gene, sep=" / ")) 
dev.off()

## Have a look by Karyotype...
C5.lm <- vector("list", length(unique(unlist(sigC5)))) # From objects to hold the model co-efficients & the r-squared values
names(C5.lm) <- unique(unlist(sigC5))
for (gn in unique(unlist(sigC5))) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat2) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~. + 0 + C5:SimpKaryo.A + C5:SimpKaryo.C) # Set the initial fully parameterised model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  C5.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
C5.out <- array(0, 
                dim=c(length(C5.lm), 3, 5), 
                dimnames=list(names(C5.lm),
                              c("C5", "SimpKaryo.A:C5", "SimpKaryo.C:C5"),
                              c("Estimate", "Std.Error", "t", "P.Value","adj.P.Value")))
for (i in 1:3) {C5.out[, i, 1:4] <- t(sapply(C5.lm, FUN=extractCoefs, coef=dimnames(C5.out)[[2]][i]))}
C5.out[,,"adj.P.Value"] <- p.adjust(C5.out[,,"P.Value"], method="fdr")
for (i in 1:3) {C5.out[which(C5.out[,i,"adj.P.Value"]>0.05),i,"Estimate"] <- 0}

C5.logFC.est <- cbind(C5.out[,"C5","Estimate"], # logFC for the heterogeneous population (M0, M3, M6, M7, Unknown)
                      C5.out[,"C5","Estimate"] + C5.out[,"SimpKaryo.A:C5","Estimate"], # logFC for M1 C5 Vs not C5
                      C5.out[,"C5","Estimate"] + C5.out[,"SimpKaryo.C:C5","Estimate"]) # logFC for M2 C5 Vs not C5
colnames(C5.logFC.est) <- c("Normal", "Abnormal", "Complex")
C5.logFC.est <- C5.logFC.est[which(rowMaxs(abs(C5.logFC.est))>0.5),] # Get rid of any where the maximum FC < sqrt(2) in all comparisons
C5.logFC.est <- C5.logFC.est[hclust(dist(C5.logFC.est))$order,] # Cluster using hclust
Gene <- sapply(rownames(C5.logFC.est), FUN=getGName, map=data.frame(ProbeID=prbData$ProbeID, TargetID=as.character(prbData$TargetID), stringsAsFactors=FALSE))
C5.ggData <- data.frame(Probe = factor(rownames(C5.logFC.est), levels=rownames(C5.logFC.est)),
                        Gene = factor(Gene, levels=unique(Gene)),
                        logFC = as.vector(C5.logFC.est),
                        SimpleKaryotype = factor(rep(colnames(C5.logFC.est), each=nrow(C5.logFC.est)), levels=c("Normal", "Abnormal","Complex")))
C5.ggData$Probe <- factor(C5.ggData$Probe, levels=C5.ggData$Probe) # Ensure the order is not changed by ggplot
C5.ggData$Gene <- factor(C5.ggData$Gene, levels=C5.ggData$Gene) # Ensure the order is not changed by ggplot
png("C5byFAB.png", height=600, width=900)
ggplot(C5.ggData, aes(x=SimpleKaryotype, y=Probe, fill=logFC)) +
  ggtitle("Estimated log2 fold-change for C5 Vs WT")+ 
  geom_raster() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(labels=paste(rownames(C5.logFC.est), Gene, sep=" / ")) 
dev.off()


## Repeat for the Core.Anchor.ID probes
coefFCA <- grep("FCA", colnames(modMat))
sigFANC.FCA <- vector("list", length(coefFCA)) #
names(sigFANC.FCA) <- colnames(modMat)[coefFCA]
for (i in coefFCA) {
  coef <- colnames(modMat)[i]
  topTab <- topTable(fit, coef=coef, number=nG)
  sigFANC.FCA[[coef]] <- rownames(topTab)[which(topTab$adj.P.Val<0.1)]
}
length(unique(unlist(sigFANC.FCA))) # [1] 172

## C1
coefC1 <- grep("C1", colnames(modMat))
sigC1 <- vector("list", length(coefC1)) #
names(sigC1) <- colnames(modMat)[coefC1]
for (i in coefC1) {
  coef <- colnames(modMat)[i]
  topTab <- topTable(fit, coef=coef, number=nG)
  sigC1[[coef]] <- rownames(topTab)[which(topTab$adj.P.Val<0.1)]
}
length(unique(unlist(sigC1))) # [1] 16

## C5
coefC5 <- grep("C5", colnames(modMat))
sigC5 <- vector("list", length(coefC5)) #
names(sigC5) <- colnames(modMat)[coefC5]
for (i in coefC5) {
  coef <- colnames(modMat)[i]
  topTab <- topTable(fit, coef=coef, number=nG)
  sigC5[[coef]] <- rownames(topTab)[which(topTab$adj.P.Val<0.1)]
}
length(unique(unlist(sigC5))) # [1] 40

