# Load libraries
library(WGCNA)

# Get info about promoter probes, defined as TSS200 or UTR5' probes
cat("Reading probes annotation file\n")
probeAnnotations <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"")
idx <- grep("TSS200|5'UTR", probeAnnotations$UCSC_REFGENE_GROUP)
promProbes <- probeAnnotations$TargetID[idx]

# Read cancer tissue methylation files (regex matches only cancer files)
cat("Retrieving file list\n")
dataFileDir = "../rawdata/GBM/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
cancerFileNames <- list.files(dataFileDir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")

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
ggCancer <- goodGenes(t(cancerBetaValues), verbose=3)

cancerMethyl <- as.data.frame(t(cancerBetaValues)[,ggCancer])
probeList <- colnames(cancerMethyl)
cancerMethylSampList <- rownames(cancerMethyl)

# Save to files
save(cancerMethyl, file="../Rdata/GBM/data/GBM-CMP.Rdata")
save(probeList,    file="../Rdata/GBM/info/probeList.Rdata")
save(cancerMethylSampList, file="../Rdata/GBM/info/methylSampList.Rdata")
quit(save="no")
