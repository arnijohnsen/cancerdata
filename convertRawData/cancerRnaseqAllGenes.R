# Set wd and load libraries
setwd("~/cancerdata/")
library(WGCNA)

# Read normal tissue rna sequencing files (regex matches only normal files)
cat("Retrieving file list\n")
dataFileDir = "../rawdata/RNASeq/RNASeq/UNC__IlluminaHiSeq_RNASeq/Level_3/"
fileNames <- list.files(dataFileDir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
# Get only annotated gene files
fileNames <- fileNames[grep("gene", fileNames)]

n <- length(fileNames)
# Read first file to get gene names and dimension of data frame
cat("Reading beta values from file", 1, "of", n, "\n")
tmpRnaseq <- read.table(paste(dataFileDir, fileNames[1], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)

genes <- gsub("\\|.*", "", tmpRnaseq$gene[-(1:29)])
# Fix one wrong names
genes[16272] <- "SLC35E2B"

rpkmValues <- data.frame(x = tmpRnaseq$RPKM[-(1:29)])
colnames(rpkmValues) <- substring(fileNames[1], 14, 27)
rownames(rpkmValues) <- genes

for (i in 2:n){
  cat("Reading beta values from file", i, "of", n, "\n")
  tmpRnaseq <- read.table(paste(dataFileDir, fileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  rpkmValues[[substring(fileNames[i], 14, 27)]] <- tmpRnaseq$RPKM[-(1:29)]
}

# Transpose and filter out bad data
gsg <- goodSamplesGenes(t(rpkmValues), verbose=3)
cancerRnaseq <- as.data.frame(t(rpkmValues)[gsg$goodSamples, gsg$goodGenes])

save(cancerRnaseq, file="../Rdata/normalRnaseqAllGenes.Rdata")
