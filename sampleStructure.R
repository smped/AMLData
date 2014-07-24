## Find the structure of the samples with respect to the other groupings
## We are really interested in the groupings as the relate to FANC & C1/C5

##########
## FANC ##
##########
struct <- list(FANC=c(), C=c())
ffnRow <- which(annot$FANC.Full.Network == "Y")
fcaRow <- which(annot$FANC.Core.Anchor.ID == "Y")
wtRow <- setdiff(1:nrow(annot), union(ffnRow, fcaRow))
struct[["FANC"]] <- matrix(0, ncol=3, nrow=42, dimnames=list(c("Source.Adelaide", "Source.ALLG", paste("M", 0:7, sep=""), 
                                                               paste("SimpKaryo", c("A","C","N"), sep="."), "PrimaryAML.Y", "PrimaryAML.N",
                                                               "FLT3.ITD", "FLT3.TKD", "NPM1", "CEBPA", "GADD45.Low", "GADD45.High", "KLF5.Low",
                                                               "KLF5.High", "Seq.DNMT3a", "Seq.FLT3.TKD", "Seq.IDH1", "Seq.IDH2", "Seq.NPM1",
                                                               "Seq.NRAS", "Seq.WT1", "FANC.Full.Network", "FANC.Core.Anchor.ID", "FANC.ALL", 
                                                               "MBD.Full.Network", "MBD.s.and.Likes", "MBD1.4.MECP2", "METC", paste("C",1:5, sep="")),
                                                             c("FANC.Full.Network", "FANC.Core.Anchor.ID", "WT")
                                                             ))
## Populate the table
tRow = grep("Adelai", annot$Source)
struct$FANC["Source.Adelaide",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("ALLG", annot$Source)
struct$FANC["Source.ALLG",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M0", annot$Phenotype)
struct$FANC["M0",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M1", annot$Phenotype)
struct$FANC["M1",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M2", annot$Phenotype)
struct$FANC["M2",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M3", annot$Phenotype)
struct$FANC["M3",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M4", annot$Phenotype)
struct$FANC["M4",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M5", annot$Phenotype)
struct$FANC["M5",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M6", annot$Phenotype)
struct$FANC["M6",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("M7", annot$Phenotype)
struct$FANC["M7",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("A", annot$SimpKaryo)
struct$FANC["SimpKaryo.A",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("C", annot$SimpKaryo)
struct$FANC["SimpKaryo.C",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("N", annot$SimpKaryo)
struct$FANC["SimpKaryo.N",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$PrimaryAML)
struct$FANC["PrimaryAML.Y",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("N", annot$PrimaryAML)
struct$FANC["PrimaryAML.N",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$FLT3.ITD)
struct$FANC["FLT3.ITD",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$FLT3.TKD)
struct$FANC["FLT3.TKD",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$NPM1)
struct$FANC["NPM1",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("b", annot$CEBPA)
struct$FANC["CEBPA",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("low", annot$GADD45..Simple.)
struct$FANC["GADD45.Low",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("high", annot$GADD45..Simple.)
struct$FANC["GADD45.High",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("LOW", annot$KLF5..Simple.)
struct$FANC["KLF5.Low",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("HIGH", annot$KLF5..Simple.)
struct$FANC["KLF5.High",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Mutant", annot$Seq.DNMT3a)
struct$FANC["Seq.DNMT3a",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Mutant", annot$Seq.FLT3.TKD)
struct$FANC["Seq.FLT3.TKD",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Mutant", annot$Seq.IDH1)
struct$FANC["Seq.IDH1",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Mutant", annot$Seq.IDH2)
struct$FANC["Seq.IDH2",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Mutant", annot$Seq.NPM1)
struct$FANC["Seq.NPM1",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Mutant", annot$Seq.NRAS)
struct$FANC["Seq.NRAS",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Mutant", annot$Seq.WT1)
struct$FANC["Seq.WT1",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$FANC.Full.Network)
struct$FANC["FANC.Full.Network",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$FANC.Core.Anchor.ID)
struct$FANC["FANC.Core.Anchor.ID",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$FANC.ALL)
struct$FANC["FANC.ALL",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$MBD.Full.Network)
struct$FANC["MBD.Full.Network",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$MBD.s.and.Likes)
struct$FANC["MBD.s.and.Likes",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$MBD1.4.MECP2)
struct$FANC["MBD1.4.MECP2",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$METC)
struct$FANC["METC",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$C1)
struct$FANC["C1",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$C2)
struct$FANC["C2",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$C3)
struct$FANC["C3",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$C4)
struct$FANC["C4",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
tRow = grep("Y", annot$C5)
struct$FANC["C5",] = c(length(intersect(tRow, ffnRow)), length(intersect(tRow, fcaRow)), length(intersect(tRow, wtRow)))
struct$FANC
## Points of note:
## 1 - All FCA mutations also have FFN mutations
## 2 - To be ignored for all analyses: M6, M7, FLT3.TKD, CEBPA, Seq.FLT3.TKD, Seq.NRAS, C2-C4
## 3 - For FCA, the annotations M0 & M3-M7 contain too few samples, however when looking at FFN alone, M0 & M4-M5 may be useful
## 4 - Potential interactions based purely on there being large enough samples are:
##                  FFN:M0, FFN:M1, FFN:M2, FFN:M4, FFN:M5, FFN:SimpKaryo(A/C/N), FFN:C1

## Another way to look would be to check the % of each annotation which are also FFN as this may give some idea about any aliasing:

###################################################################################################################
###################################################################################################################
##                                       Just do a correlation plot                                              ##
###################################################################################################################
###################################################################################################################

## The vectors/objects we can already use are Source, M, SimpKaryo, PrimaryAML, NPM1, GADD45, FANC, C1, C5
## The rest need to be formed
FLT3 <- matrix(0, nrow=n, ncol=2)
colnames(FLT3) <- c("FLT3.ITD", "FLT3.TKD")
FLT3[which(annot$FLT3.ITD=="Y"),"FLT3.ITD"] <- 1
FLT3[which(annot$FLT3.TKD=="Y"),"FLT3.TKD"] <- 1

KLF5 <- matrix(0, nrow=n, ncol=2)
colnames(KLF5) <- c("KLF5.Low", "KLF5.High")
KLF5[grep("HIGH",annot$KLF5..Simple.), "KLF5.High"] <- 1
KLF5[grep("LOW",annot$KLF5..Simple.), "KLF5.Low"] <- 1

Seq <- matrix(0, nrow=n, ncol=7)
colnames(Seq) <- paste("Seq", c("DNMT3a", "FLT3.TKD", "IDH1", "IDH2", "NPM1", "NRAS", "WT1"), sep=".")
for (i in 1:7) Seq[which(annot[,colnames(Seq)[i]]=="Mutant"), i] <- 1

MBD <- matrix(0, nrow=n, ncol=3)
colnames(MBD) <- paste("MBD", c(".Full.Network", ".s.and.Likes", "1.4.MECP2"), sep="")
for (i in 1:3) MBD[which(annot[,colnames(MBD)[i]]=="Y"),i] <- 1

METC <- rep(0, n)
METC[which(annot$METC=="Y")] <- 1

fullMat <- cbind(Source, M, SimpKaryo, PrimaryAML, FLT3, NPM1, GADD45, KLF5, Seq, FANC[,1:2], MBD, METC, C1, C5)
mcor <- cor(fullMat)
png("correlations.png", height=800, width=800)
corrplot(mcor)
dev.off()
## The good news is that there should be minimal aliasing with FANC & any other predictors!!!

## Calculate the weights based on the full matrix
arrWeights <- arrayWeights(exprs(normPrbData)[detPrbs,useArr], design=fullMat, weights=prbWeights[detPrbs,useArr])
