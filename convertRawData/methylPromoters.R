# Set wd and load libraries
setwd("~/cancerdata/")
library(WGCNA)

# Get info about promoter probes, defined as TSS200 or UTR5' probes
cat("Reading probes annotation file\n")
probeAnnotations <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"")
idx <- grep("TSS200|5'UTR", probeAnnotations$UCSC_REFGENE_GROUP)
promProbes <- probeAnnotations$TargetID[idx]

# Read normal tissue methylation files (regex matches only normal files)
cat("Retrieving file list\n")
dataFileDir = "../rawdata/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
normalFileNames <- list.files(dataFileDir, pattern=".*1[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
cancerFileNames <- list.files(dataFileDir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")

# Create new data.frame containing beta values, 
# with each normal sample in one column and each probe in one row
normalBetaValues <- data.frame(rownames = promProbes)
nNormal <- length(normalFileNames)
cat("Reading normal beta values\n")
pb <- txtProgressBar(min=1, max=nNormal, style=3)
for (i in 1:nNormal){
  setTxtProgressBar(pb, i)
  tmpMethyl <- read.table(paste(dataFileDir, normalFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmpMethyl[tmpMethyl$Composite.Element.REF %in% promProbes, ]
  normalBetaValues[[substring(gsub(".*TCGA", "TCGA", normalFileNames[i]),1,14)]] <- tmp$Beta_value 
}
cat("\n")
normalBetaValues$rownames <- NULL
rownames(normalBetaValues) <- promProbes

# Create new data.frame containing beta values, 
# with each cancer sample in one column and each probe in one row
cancerBetaValues <- data.frame(rownames = promProbes)
nCancer <- length(cancerFileNames)
cat("Reading cancer beta values\n")
pb <- txtProgressBar(min=1, max=nCancer, style=3)
for (i in 1:nCancer){
  setTxtProgressBar(pb, i)
  tmpMethyl <- read.table(paste(dataFileDir, cancerFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmpMethyl[tmpMethyl$Composite.Element.REF %in% promProbes, ]
  cancerBetaValues[[substring(gsub(".*TCGA", "TCGA", cancerFileNames[i]),1,14)]] <- tmp$Beta_value 
}
cat("\n")
cancerBetaValues$rownames <- NULL
rownames(cancerBetaValues) <- promProbes

# Transpose and filter out bad data
ggNormal <- goodGenes(t(normalBetaValues), verbose=3)
ggCancer <- goodGenes(t(cancerBetaValues), verbose=3)

normalMethyl <- as.data.frame(t(normalBetaValues)[, ggNormal & ggCancer])
cancerMethyl <- as.data.frame(t(cancerBetaValues)[, ggNormal & ggCancer])
probeList <- colnames(normalMethyl)
normalMethylSampList <- rownames(normalMethyl)
cancerMethylSampList <- rownames(cancerMethyl)

# Save to files
save(normalMethyl, file="../Rdata/normalMethylPromoters.Rdata")
save(cancerMethyl, file="../Rdata/cancerMethylPromoters.Rdata")
save(probeList,    file="../Rdata/probeList.Rdata")
save(normalMethylSampList, cancerMethylSampList, file="../Rdata/methylSampList.Rdata")
