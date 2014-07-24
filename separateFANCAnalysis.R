source("loadPackages.R")
source("extraFunctions.R")

## Form a model matrix:
modMat <- cbind(1, Source, M, PrimaryAML, FLT3, NPM1, GADD45, FANC[,1:2], MBD, C1, C5)
colnames(modMat)[1] <- "Intercept"
rownames(modMat) <- annot[,"MicroarrayIDinDB"]

## Estimate the weights 
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

## Check the FANC.Full.Network probes
topTab.FFN <- topTable(fit, coef="FANC.Full.Network", adjust.method="fdr", number=nG, sort.by="P")
length(which(topTab.FFN$adj.P.Val < 0.05)) # [1] 82
length(which(abs(topTab.FFN$logFC[1:82])>0.5)) # [1] 39
## OK, that actually looks better than when they were merged

## Inspect the Core.Anchor.ID probes
topTab.FCA <- topTable(fit, coef="FANC.Core.Anchor.ID", adjust.method="fdr", number=nG, sort.by="P")
length(which(topTab.FCA$adj.P.Val < 0.05)) # [1] 304
length(which(abs(topTab.FCA$logFC[1:304])>0.5)) # [1] 148
## Looks like this was a far better approach

## Try to avoid overfitting again
## FANC.Full.Network
sigFANC.FFN <- rownames(topTab.FFN)[which(abs(topTab.FFN$logFC)[1:82]>0.5)]
modCols <- apply(fit$p.value[sigFANC.FFN,], FUN=function(x,names){names[which(x<0.01)]}, MARGIN=1, names=colnames(modMat))

n <- length(useArr)
cookPoints.FFN <- vector("list", length(sigFANC.FFN)); names(cookPoints.FFN) <- sigFANC.FFN # Holds any outliers using Cook's distance
lmFANC.FFNpre <- vector("list", length(sigFANC.FFN)); names(lmFANC.FFNpre) <- sigFANC.FFN # Holds the model coefficients before outlier omission
lmFANC.FFNpost <- lmFANC.FFNpre
r2FANC.FFN <- vector("list", length(sigFANC.FFN)); names(r2FANC.FFN) <- sigFANC.FFN # Holds the R^2 values
for (gn in sigFANC.FFN) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat[,modCols[[gn]]])
  fm <- update(formula(tempDat),~.+0)
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr])
  lmFANC.FFNpre[[gn]] <- summary(tempLm)$coef
  p <- length(modCols[[gn]]) # The number of parameters being fitted
  cp <- which(cooks.distance(tempLm) > qf(0.1, p, n-p )) # Find any overly influential points
  if (length(cp) > 0) {
    tempLm <- lm(fm, data=tempDat[-cp,], weights=wtMat[gn, useArr[-cp]])
    cookPoints.FFN[[gn]] <- cp
  }
  lmFANC.FFNpost[[gn]] <- summary(tempLm)$coef
  r2FANC.FFN[[gn]] <- summary(tempLm)$r.squared
}
## Check the models after omission of outliers
postFANC.FFN <- t(data.frame(lapply(lmFANC.FFNpost, FUN = function(x) x["FANC.Full.Network",])))
postFANC.FFN <- data.frame(Probe=substr(rownames(postFANC.FFN), start=2, stop=nchar(rownames(postFANC.FFN))), postFANC.FFN, adj.P=p.adjust(postFANC.FFN[,4], method="fdr"))
postFANC.FFN %.% arrange(adj.P) %.% filter(adj.P<0.05)
## And before the omission of outliers
preFANC.FFN <- t(data.frame(lapply(lmFANC.FFNpre, FUN = function(x) x["FANC.Full.Network",])))
preFANC.FFN <- data.frame(Probe=substr(rownames(postFANC.FFN), start=2, stop=nchar(rownames(postFANC.FFN))), preFANC.FFN, adj.P=p.adjust(preFANC.FFN[,4], method="fdr"))
preFANC.FFN %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)
## Maybe just use the pre-removal dataset & just note the overly influential points


sigFANC.FCA <- rownames(topTab.FCA)[which(abs(topTab.FCA$logFC)[1:467]>0.5)]
modCols <- apply(fit$p.value[sigFANC.FCA,], FUN=function(x,names){names[which(x<0.01)]}, MARGIN=1, names=colnames(modMat))

n <- length(useArr)
cookPoints.FCA <- vector("list", length(sigFANC.FCA)); names(cookPoints.FCA) <- sigFANC.FCA # Holds any outliers using Cook's distance
lmFANC.FCApre <- vector("list", length(sigFANC.FCA)); names(lmFANC.FCApre) <- sigFANC.FCA # Holds the model coefficients before outlier omission
lmFANC.FCApost <- lmFANC.FCApre
r2FANC.FCA <- vector("list", length(sigFANC.FCA)); names(r2FANC.FCA) <- sigFANC.FCA # Holds the R^2 values
for (gn in sigFANC.FCA) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat[,modCols[[gn]]])
  fm <- update(formula(tempDat),~.+0)
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr])
  lmFANC.FCApre[[gn]] <- summary(tempLm)$coef
  p <- length(modCols[[gn]]) # The number of parameters being fitted
  cp <- which(cooks.distance(tempLm) > qf(0.1, p, n-p )) # Find any overly influential points
  if (length(cp) > 0) {
    tempLm <- lm(fm, data=tempDat[-cp,], weights=wtMat[gn, useArr[-cp]])
    cookPoints.FCA[[gn]] <- cp
  }
  lmFANC.FCApost[[gn]] <- summary(tempLm)$coef
  r2FANC.FCA[[gn]] <- summary(tempLm)$r.squared
}
## Check the models after omission of outliers
postFANC.FCA <- t(data.frame(lapply(lmFANC.FCApost, FUN = function(x) x["FANC.Core.Anchor.ID",])))
postFANC.FCA <- data.frame(Probe=substr(rownames(postFANC.FCA), start=2, stop=nchar(rownames(postFANC.FCA))), postFANC.FCA, adj.P=p.adjust(postFANC.FCA[,4], method="fdr"))
postFANC.FCA %.% arrange(adj.P) %.% filter(adj.P<0.05)
## And before the omission of outliers
preFANC.FCA <- t(data.frame(lapply(lmFANC.FCApre, FUN = function(x) x["FANC.Core.Anchor.ID",])))
preFANC.FCA <- data.frame(Probe=substr(rownames(postFANC.FCA), start=2, stop=nchar(rownames(postFANC.FCA))), preFANC.FCA, adj.P=p.adjust(preFANC.FCA[,4], method="fdr"))
preFANC.FCA %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)
## Maybe just use the pre-removal dataset & just note the overly influential points

## Overall, the models before removal of influential points seem to be more sensible. 
## It's impossible to tell if those points are good or bad data points, so leave them in & just let the data speak.

## Use the preFANC.FCA & FFN datasets. Map the probes to genes & send them to Mahmoud.
preFANC.FCA <- data.frame(Probe=preFANC.FCA$Probe, Gene=sapply(as.character(preFANC.FCA$Probe),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}),
                          preFANC.FCA[,2:4], rawP = preFANC.FCA[,5], adj.P = preFANC.FCA$adj.P)
preFANC.FCA %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)
write.csv(preFANC.FCA %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1), "FANC.Core.Anchor.ID.Probes.csv", row.names=FALSE)

## Use the preFANC.FFN & FFN datasets. Map the probes to genes & send them to Mahmoud.
preFANC.FFN <- data.frame(Probe=preFANC.FFN$Probe, Gene=sapply(as.character(preFANC.FFN$Probe),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}),
                          preFANC.FFN[,2:4], rawP = preFANC.FFN[,5], adj.P = preFANC.FFN$adj.P)
preFANC.FFN %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)
write.csv(preFANC.FFN %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1), "FANC.Full.Network.Probes.csv", row.names=FALSE)

#############
## C1 & C5 ##
#############

## Inspect the C1
topTab.C1 <- topTable(fit, coef="C1", adjust.method="fdr", number=nG, sort.by="P")
length(which(topTab.C1$adj.P.Val < 0.05)) # [1] 92
length(which(abs(topTab.C1$logFC[1:92])>0.5)) # [1] 28

sigC1 <- rownames(topTab.C1)[which(abs(topTab.C1$logFC)[1:92]>0.5)]
modCols <- apply(fit$p.value[sigC1,], FUN=function(x,names){names[which(x<0.01)]}, MARGIN=1, names=colnames(modMat))

lmC1 <- vector("list", length(sigC1)); names(lmC1) <- sigC1 # Holds the model coefficients before outlier omission
r2C1 <- vector("list", length(sigC1)); names(r2C1) <- sigC1 # Holds the R^2 values
for (gn in sigC1) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat[,modCols[[gn]]])
  fm <- update(formula(tempDat),~.+0)
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr])
  rTerms <- names(which(summary(tempLm))
  lmC1[[gn]] <- summary(tempLm)$coef
  r2C1[[gn]] <- summary(tempLm)$r.squared
}
## And before the omission of outliers
C1.out <- t(data.frame(lapply(lmC1, FUN = function(x) x["C1",])))
C1.out <- data.frame(Probe=substr(rownames(C1.out), start=2, stop=nchar(rownames(C1.out))), C1.out, adj.P=p.adjust(C1.out[,4], method="fdr"))
C1.out %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)

## Use the C1.out datasets. Map the probes to genes & send them to Mahmoud.
C1.out <- data.frame(Probe=C1.out$Probe, Gene=sapply(as.character(C1.out$Probe),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}),
                     C1.out[,2:4], rawP = C1.out[,5], adj.P = C1.out$adj.P)
C1.out %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)
write.csv(C1.out %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1), "C1.Probes.csv", row.names=FALSE)

## Inspect the C5
topTab.C5 <- topTable(fit, coef="C5", adjust.method="fdr", number=nG, sort.by="P")
length(which(topTab.C5$adj.P.Val < 0.05)) # [1] 142
length(which(abs(topTab.C5$logFC[1:142])>0.5)) # [1] 53

sigC5 <- rownames(topTab.C5)[which(abs(topTab.C5$logFC)[1:142]>0.5)]
modCols <- apply(fit$p.value[sigC5,], FUN=function(x,names){names[which(x<0.01)]}, MARGIN=1, names=colnames(modMat))

lmC5 <- vector("list", length(sigC5)); names(lmC5) <- sigC5 # Holds the model coefficients before outlier omission
r2C5 <- vector("list", length(sigC5)); names(r2C5) <- sigC5 # Holds the R^2 values
for (gn in sigC5) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat[,modCols[[gn]]])
  fm <- update(formula(tempDat),~.+0)
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr])
  lmC5[[gn]] <- summary(tempLm)$coef
  p <- length(modCols[[gn]]) # The number of parameters being fitted
  r2C5[[gn]] <- summary(tempLm)$r.squared
}
## And before the omission of outliers
C5.out <- t(data.frame(lapply(lmC5, FUN = function(x) x["C5",])))
C5.out <- data.frame(Probe=substr(rownames(C5.out), start=2, stop=nchar(rownames(C5.out))), C5.out, adj.P=p.adjust(C5.out[,4], method="fdr"))
C5.out %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)



## Use the C5.out & FFN datasets. Map the probes to genes & send them to Mahmoud.
C5.out <- data.frame(Probe=C5.out$Probe, Gene=sapply(as.character(C5.out$Probe),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}),
                     C5.out[,2:4], rawP = C5.out[,5], adj.P = C5.out$adj.P)
C5.out %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1)
write.csv(C5.out %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate)>1), "C5.Probes.csv", row.names=FALSE)

###################################################################################################################################
## Try a better updating method. Under this approach the full model is fit & terms removed starting with the least significant.  ##
## Term removal stops either when all Holm's adjusted p-values are < 0.05, or the term of interest has been removed              ##
###################################################################################################################################
#######################
## FANC.Full.Network ##
#######################
sigFANC.FFN <- rownames(topTab.FFN)[which(topTab.FFN$adj.P.Val < 0.05)]
FANC.FFN.lm <- vector("list", length(sigFANC.FFN)) # From objects to hold the model co-efficients & the r-squared values
names(FANC.FFN.lm) <- sigFANC.FFN
for (gn in sigFANC.FFN) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~.+0) # Set the initial saturated model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  FANC.FFN.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
## Extract the term of interest:
FANC.FFN <- t(data.frame(lapply(FANC.FFN.lm, FUN = function(x) x["FANC.Full.Network",])))
FANC.FFN <- data.frame(Probe=sigFANC.FFN,
                       Gene=sapply(as.character(sigFANC.FFN),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}), # Get the gene from the original object
                       FANC.FFN[,1:4], 
                       adj.P=p.adjust(FANC.FFN[,4], method="holm"))
## Export only those with fold-change > 1.5
write.csv(FANC.FFN %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate) > log2(1.5)), "FANC.Full.Network.Probes.csv", row.names=FALSE)

#########################
## FANC.Core.Anchor.ID ##
#########################
sigFANC.FCA <- rownames(topTab.FCA)[which(topTab.FCA$adj.P.Val < 0.05)]
FANC.FCA.lm <- vector("list", length(sigFANC.FCA)) # Form objects to hold the model co-efficients
names(FANC.FCA.lm) <- sigFANC.FCA
for (gn in sigFANC.FCA) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~.+0) # Set the initial saturated model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  FANC.FCA.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
## Extract the term of interest:
FANC.FCA <- t(data.frame(lapply(FANC.FCA.lm, FUN = function(x) {x["FANC.Core.Anchor.ID",]})))
FANC.FCA <- data.frame(Probe=sigFANC.FCA,
                       Gene=sapply(as.character(sigFANC.FCA),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}), # Get the gene from the original object
                       FANC.FCA[,1:4], 
                       adj.P=p.adjust(FANC.FCA[,4], method="holm"))
## Export only those with fold-change > 1.5
write.csv(FANC.FCA %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate) > log2(1.5)), "FANC.Full.Network.Probes.csv", row.names=FALSE)

########
## C1 ##
########
sigC1 <- rownames(topTab.C1)[which(topTab.C1$adj.P.Val < 0.05)]
C1.lm <- vector("list", length(sigC1)) # Form objects to hold the model co-efficients
names(C1.lm) <- sigC1
for (gn in sigC1) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~.+0) # Set the initial saturated model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  C1.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
## Extract the term of interest:
C1 <- t(data.frame(lapply(C1.lm, FUN = function(x) {x["C1",]})))
C1 <- data.frame(Probe=sigC1,
                       Gene=sapply(as.character(sigC1),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}), # Get the gene from the original object
                       C1[,1:4], 
                       adj.P=p.adjust(C1[,4], method="holm"))
## Export only those with fold-change > 1.5
write.csv(C1 %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate) > log2(1.5)), "C1.Probes.csv", row.names=FALSE)

########
## C5 ##
########
sigC5 <- rownames(topTab.C5)[which(topTab.C5$adj.P.Val < 0.05)]
C5.lm <- vector("list", length(sigC5)) # Form objects to hold the model co-efficients
names(C5.lm) <- sigC5
for (gn in sigC5) {
  tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat) # Form a temporary data object with the probe data & the full model matrix
  fm <- update(formula(tempDat),~.+0) # Set the initial saturated model
  tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Evaluate the initial model
  C5.lm[[gn]] <- summary(step(tempLm))$coef # Keep the final model using an automated step-wise model selection process
}
## Extract the term of interest:
C5 <- t(data.frame(lapply(C5.lm, FUN = function(x) {x["C5",]})))
C5 <- data.frame(Probe=sigC5,
                 Gene=sapply(as.character(sigC5),FUN=function(x){as.character(prbData[which(prbData[,2]==x),1])}), # Get the gene from the original object
                 C5[,1:4], 
                 adj.P=p.adjust(C5[,4], method="holm"))
## Export only those with fold-change > 1.5
write.csv(C5 %.% arrange(adj.P) %.% filter(adj.P<0.05, abs(Estimate) > log2(1.5)), "C5.Probes.csv", row.names=FALSE)

###########################################################
###########################################################
##              Heatmap plotting                         ##
###########################################################
###########################################################

## Group the arrays in a specific order with the C5 samples at the far left
o <- order(modMat[,"C5"], modMat[,"M5"], modMat[,"M4"], modMat[,"M0"], modMat[,"FANC.Core.Anchor.ID"], modMat[,"Source"])

## Plot a heatmap for each of the 4 sets of genes declared significant
C5sigProbes <- as.character(C5$Probe[intersect(which(C5$adj.P < 0.05), which(abs(C5$Estimate)>log2(1.5))), drop=TRUE])
## Hmmm. This is tricky...

######################################
## Maybe check interactions as well ##
######################################

## e.g 
tempDat <- data.frame(exprs = exprs(normPrbData)[gn, useArr], modMat)
fm <- update(formula(tempDat),~(.+0+M0:C5+M1:C5+M2:C5+M3:C5+M4:C5+M5:C5)) ## Incorporate the various FAB stages interacting with C5 
tempLm <- lm(fm, data=tempDat, weights=wtMat[gn, useArr]) # Eval
finLm <- step(tempLm)
summary(finLm)
