## A very simple look at the data just comparing FANC with non-FANC as a two way comparison

modMatFFNvsWT <- cbind(FANC[,1], 1-FANC[,1])
colnames(modMatFFNvsWT) <- c("FFN", "WT")
rownames(modMatFFNvsWT) <- useArr

## Estimate the weights without the probe weights
arrWeightsFFNvsWT <- arrayWeights(exprs(normPrbData)[detPrbs,useArr], design=modMatFFNvsWT)#, weights=prbWeights[detPrbs,useArr])
barplot(arrWeightsFFNvsWT, names.arg=useArr, las=2)
abline(h=1, col="blue", lty=2)
wtMatFFNvsWT <- arrWeightsFFNvsWT# * prbWeights[detPrbs,useArr]

## Fit the model
fitFFNvsWT <- lmFit(exprs(normPrbData)[detPrbs,useArr], design=modMatFFNvsWT, weights=wtMatFFNvsWT)
contMatFFNvsWT <- makeContrasts(FFN-WT, levels=modMatFFNvsWT)
fitFFNvsWT <- contrasts.fit(fitFFNvsWT, contMatFFNvsWT)
fitFFNvsWT <- eBayes(fitFFNvsWT)
## get the top FANC genes
topTabFFNvsWT <- topTable(fitFFNvsWT, n=nG)
length(which(topTabFFNvsWT$adj.P.Val < 0.05)) # 14
sigFFNvsWT <- rownames(topTabFFNvsWT)[intersect(which(topTabFFNvsWT$adj.P.Val < 0.05), which(abs(topTabFFNvsWT$logFC) > 1))]
length(sigFFNvsWT) # [1] 0
## It's the weights that give the false readings!!!


## Manual inspection is showing all of these to be bullshit!!!
curG <- sigFFNvsWT[1] # Pick a number they're all crap...
lapply(list(WT=exprs(normPrbData)[curG, useArr[which(FANC[,1]==0)]], FANC=exprs(normPrbData)[curG, useArr[which(FANC[,1]==1)]]), FUN=mean)
boxplot(list(WT=exprs(normPrbData)[curG, useArr[which(FANC[,1]==0)]], FANC=exprs(normPrbData)[curG, useArr[which(FANC[,1]==1)]]))
## I have NFI what is going wrong. All I can think of is there is one or two massively weighted points which are skewing things

