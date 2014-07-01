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
normFileNames <- list.files(dataFileDir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")

# Create new data.frame containing beta values, 
# with each normal sample in one column and each probe in one row
betaValues <- data.frame(rownames = promProbes)
n <- length(normFileNames)
for (i in 1:n){
  cat("Reading beta values from file", i, "of", n, "\n")
  tmpMethyl <- read.table(paste(dataFileDir, normFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmpMethyl[tmpMethyl$Composite.Element.REF %in% promProbes, ]
  betaValues[[substring(normFileNames[i],47,60)]] <- tmp$Beta_value 
}
betaValues$rownames <- NULL
rownames(betaValues) <- promProbes

# Transpose and filter out bad data
gsg <- goodSamplesGenes(t(betaValues), verbose=3)
cancerMethyl <- as.data.frame(t(betaValues)[gsg$goodSamples, gsg$goodGenes])

save(cancerMethyl, file="../Rdata/cancerMethylPromoters.Rdata")
