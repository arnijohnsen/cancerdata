# Load libraries
library(WGCNA)

# Get info about promoter probes, defined as TSS200 or UTR5' probes
cat("Reading probes annotation file\n")
probeAnnotations <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"")
idxZone   <- !grepl("TSS200|TSS1500|1stExon|5'UTR|3'UTR|^$", probeAnnotations$UCSC_REFGENE_GROUP)
idxIsland <-  grepl("Island", probeAnnotations$RELATION_TO_UCSC_CPG_ISLAND)
islandProbes <- probeAnnotations$TargetID[idxZone & idxIsland]

# Read normal tissue methylation files (regex matches only normal files)
cat("Retrieving file list\n")
dataFileDir = "../rawdata/BRCA/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
normalFileNames <- list.files(dataFileDir, pattern=".*1[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
cancerFileNames <- list.files(dataFileDir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")

# Create new data.frame containing beta values, 
# with each normal sample in one column and each probe in one row
normalBetaValues <- data.frame(rownames = islandProbes)
nNormal <- length(normalFileNames)
cat("Reading normal beta values\n")
pb <- txtProgressBar(min=1, max=nNormal, style=3)
for (i in 1:nNormal){
  setTxtProgressBar(pb, i)
  tmpMethyl <- read.table(paste(dataFileDir, normalFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmpMethyl[tmpMethyl$Composite.Element.REF %in% islandProbes, ]
  normalBetaValues[[substring(gsub(".*TCGA", "TCGA", normalFileNames[i]),1,14)]] <- tmp$Beta_value
}
cat("\n")
normalBetaValues$rownames <- NULL
rownames(normalBetaValues) <- islandProbes

# Create new data.frame containing beta values, 
# with each cancer sample in one column and each probe in one row
cancerBetaValues <- data.frame(rownames = islandProbes)
nCancer <- length(cancerFileNames)
cat("Reading cancer beta values\n")
pb <- txtProgressBar(min=1, max=nCancer, style=3)
for (i in 1:nCancer){
  setTxtProgressBar(pb, i)
  tmpMethyl <- read.table(paste(dataFileDir, cancerFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmpMethyl[tmpMethyl$Composite.Element.REF %in% islandProbes, ]
  cancerBetaValues[[substring(gsub(".*TCGA", "TCGA", cancerFileNames[i]),1,14)]] <- tmp$Beta_value 
}
cat("\n")
cancerBetaValues$rownames <- NULL
rownames(cancerBetaValues) <- islandProbes

# Save data for next script, as doing everything in the same script
# caused memory issues
save(normalBetaValues, file="../Rdata/BRCA/tmp/normalBetaValuesIsland.Rdata")
save(cancerBetaValues, file="../Rdata/BRCA/tmp/cancerBetaValuesIsland.Rdata")
