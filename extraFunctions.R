## Functions for combining p-values

Stouffer <- function(x) { # p is a vector of p-values
  ## x must be a matrix with p-values in the first column & weights in the second
  Zi <- qnorm(1-x[,1]) 
  Z  <- sum(x[,2]*Zi)/sqrt(sum(x[,2]^2))
  p.val <- 1-pnorm(Z)
  return(p.val)
}

## Extract the key coefficients from each probe
extractCoefs <- function (summary, coef){
  rw <- which(rownames(summary)==coef)
  if (length(rw) == 1) return(summary[rw,])
  else return(c(rep(0, 3), NA))
}

## Extract the key coefficients from each probe
extractCoefs.lm <- function (lm, coef){
  summary <- summary(lm)$coef
  rw <- which(rownames(summary)==coef)
  if (length(rw) == 1) return(summary[rw,])
  else return(c(rep(0, 3), NA))
}

aster <- function(x) {
  out <- rep("", length(x))
  out[which(x<0.1)] <- "."
  out[which(x<0.05)] <- "*"
  out[which(x<0.01)] <- "**"
  out[which(x<0.001)] <- "***"
  out
}

## Get the gene names form the inital mapping
getGName <- function(probe, map) {
  rw <- which(map[,"ProbeID"]==probe)
  map[rw,"TargetID"]
}

## Rotate a txt file...
transposeTxt <- function(inFile, outFile=NULL, ...) {
  ## Will export the file if the outFile is specified
  ## Otherwise will just load it into memory
  if (!file.exists(inFile)) return(paste(inFile, " could not be found", sep=""))
  x <- read.delim(inFile, stringsAsFactors=FALSE, header=FALSE)
  out <- data.frame(t(x), stringsAsFactors=FALSE)
  if (!is.null(outFile)) {
    write.table(out, outFile, ...)
    return(paste("File written to ", outFile, sep=""))
  }
  else out
}