source("loadPackages.R")
source("extraFunctions.R")

## Organise the annotations. This will take a while.
## NB: I've edited this file outside of R so it's in the correct orientation.
annotFile <- file.path("data", "FullAnnotation.csv")
annot <- read.csv(annotFile, header=TRUE, stringsAsFactors=FALSE)
## Correct the sample names
annot[,1] <- paste("X_", substr(paste("00", annot[,1], sep=""), start=nchar(annot[,1]), stop=nchar(annot[,1])+2), sep="")

###########################################
## Remove the Unknowns where appropriate ##
###########################################

############
## Source ##
############
annot[which(annot$Source=="Unknown"),]
## The "Unknown" source samples have a much poorer level of annotation across all columns, so remove these
## They may be useful if looking at the FLT3 status, but that's about it

#########
## Sex ##
#########
annot[which(annot$Sex==0),]

###################
## FAB groupings ##
###################
annot[which(annot$Phenotype=="Unknown"),]
## One of them is actually informative (158) which is number 6, but remove it anyway
annot$PrimaryAML[which(annot$Phenotype=="CD34")] <- "N"

######################
## Simple Karyotype ##
######################
annot[which(annot$SimpKaryo=="Unknown"),]
## Two of these contain no information about FANC, but leave them in anyway

## Just remove the most poorly annotated
rem <- unique(c(which(annot$Source=="Unknown"), which(annot$Sex==0), which(annot$SimpKaryo=="Unknown")))
annot <- annot[-rem, -7] # We don't need the detialed karyotype column either...
n <- nrow(annot) # [1] That's 143 now so we've lost 25
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

## For the PrimaryAML, we can also take any CD34 entries from the 
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

## There are only 1, 1 & 2 samples from C2, C3 & C4 so these columns are a little redundant
C1 <- C5 <- rep(0, n)
C1[which(annot$C1=="Y")] <- 1
C5[which(annot$C5=="Y")] <- 1

## Form a model matrix:
modMat <- cbind(1, Source, M, PrimaryAML, FLT3.Merged, NPM1, GADD45, FANC.Merged, MBD, C1)
colnames(modMat)[1] <- "Intercept"
rownames(modMat) <- annot[,"MicroarrayIDinDB"]

## Estimate the weights 
useArr <- samples[-rem] ## The included arrays
arrWeights <- arrayWeights(exprs(normPrbData)[detPrbs,useArr], design=modMat, weights=prbWeights[detPrbs,useArr])
barplot(arrWeights, names.arg=useArr, las=2)
abline(h=1, col="blue", lty=2)
wtMat <- arrWeights * prbWeights[detPrbs,useArr]

## Fit the model
fit <- lmFit(exprs(normPrbData)[detPrbs,useArr], design=modMat, weights=wtMat)
fit <- eBayes(fit)
## Adjust as an overall set of p-values, but excluding the intercept from the adjustment as they are all significant (by implication after DABG)
fit.adj <- fit
fit.adj$p.value[,-1] <- matrix(p.adjust(fit$p.value[,-1], method="fdr"), ncol=ncol(fit$p.value[,-1]), byrow=FALSE) 
dimnames(fit.adj$p.value) <- dimnames(fit$p.value)
topTabFANC <- topTable(fit.adj, coef="FANC.Merged", adjust.method="none", number=nG, sort.by="P")
length(which(topTabFANC$adj.P.Val<0.05)) # [1] 467
length(which(abs(topTabFANC$logFC)[1:467]>0.5)) # [1] 58
## That gives 58 probes to investigate further. Fit a model manually to avoid over-fitting
sigFANC <- rownames(topTabFANC)[which(abs(topTabFANC$logFC)[1:467]>0.5)]

## Manually remove terms from the model based on the p-values. Include terms based on a raw p-value < 0.01.
## This will avoid a degree of over-fitting but still incorporate the terms with the strongest evidence for having an effect
modCols <- apply(fit$p.value[sigFANC,], FUN=function(x,names){names[which(x<0.01)]}, MARGIN=1, names=colnames(modMat))

## We can fit models for each of the significant genes, but it is quite likely that there will be highly influential points
## The process should be:
## 1 - Fit the model with the significant terms & using the combined weights
## 2 - Find any points with Cook's distance > 1 & note them
## 3 - Fit the model with those points removed & see if the FANC terms are still significant & with |logFC| > 1

n <- length(useArr)
cookPoints <- vector("list", length(sigFANC)); names(cookPoints) <- sigFANC # Holds any outliers using Cook's distance
lmFANCpre <- vector("list", length(sigFANC)); names(lmFANCpre) <- sigFANC # Holds the model coefficients before outlier omission
lmFANCpost <- lmFANCpre
r2FANC <- vector("list", length(sigFANC)); names(r2FANC) <- sigFANC # Holds the R^2 values
for (gn in sigFANC) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat[,modCols[[gn]]])
  fm <- update(formula(tempDat),~.+0)
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr])
  lmFANCpre[[gn]] <- summary(tempLm)$coef
  p <- length(modCols[[gn]]) # The number of parameters being fitted
  cp <- which(cooks.distance(tempLm) > qf(0.1, p, n-p )) # Find any overly influential points
  if (length(cp) > 0) {
    tempLm <- lm(fm, data=tempDat[-cp,], weights=wtMat[gn, useArr[-cp]])
    cookPoints[[gn]] <- cp
  }
  lmFANCpost[[gn]] <- summary(tempLm)$coef
  r2FANC[[gn]] <- summary(tempLm)$r.squared
}
## Check the models after omission of outliers
postFANC <- t(data.frame(lapply(lmFANCpost, FUN = function(x) x["FANC.Merged",])))
postFANC <- data.frame(Probe=substr(rownames(postFANC), start=2, stop=nchar(rownames(postFANC))), postFANC, adj.P=p.adjust(postFANC[,4], method="fdr"))
## And before the omissino of outliers
preFANC <- t(data.frame(lapply(lmFANCpre, FUN = function(x) x["FANC.Merged",])))
preFANC <- data.frame(preFANC, adj.P=p.adjust(preFANC[,4], method="fdr"))
preFANC <- preFANC[order(preFANC[,4]),]

## Some of those take off in opposite directions after removal of outliers. Maybe this is the separate FANC mutations at play
## Move to the next stage which is to separate out the FANC terms...


