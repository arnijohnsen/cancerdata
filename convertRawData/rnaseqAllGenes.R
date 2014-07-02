# Set wd and load libraries
setwd("~/cancerdata/")
library(WGCNA)

# Read normal tissue rna sequencing files (regex matches only normal files)
cat("Retrieving file list\n")
dataFileDir = "../rawdata/RNASeq/RNASeq/UNC__IlluminaHiSeq_RNASeq/Level_3/"
normalFileNames <- list.files(dataFileDir, pattern=".*1[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
cancerFileNames <- list.files(dataFileDir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
# Get only annotated gene files
normalFileNames <- normalFileNames[grep("gene", normalFileNames)]
cancerFileNames <- cancerFileNames[grep("gene", cancerFileNames)]

nNormal <- length(normalFileNames)
# Read first file to get gene names and dimension of data frame
cat("Reading normal rpkm values from file", 1, "of", nNormal, "\n")
tmpRnaseq <- read.table(paste(dataFileDir, normalFileNames[1], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmpRnaseq$gene[-(1:29)])
# Fix one wrong names
genes[16272] <- "SLC35E2B"
normalRpkmValues <- data.frame(x = tmpRnaseq$RPKM[-(1:29)])
colnames(normalRpkmValues) <- substring(normalFileNames[1], 14, 27)
rownames(normalRpkmValues) <- genes
for (i in 2:nNormal){
  cat("Reading normal rpkm values from file", i, "of", nNormal, "\n")
  tmpRnaseq <- read.table(paste(dataFileDir, normalFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  normalRpkmValues[[substring(normalFileNames[i], 14, 27)]] <- tmpRnaseq$RPKM[-(1:29)]
}

nCancer <- length(cancerFileNames)
# Read first file to get gene names and dimension of data frame
cat("Reading cancer rpkm values from file", 1, "of", nCancer, "\n")
tmpRnaseq <- read.table(paste(dataFileDir, cancerFileNames[1], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmpRnaseq$gene[-(1:29)])
# Fix one wrong names
genes[16272] <- "SLC35E2B"
cancerRpkmValues <- data.frame(x = tmpRnaseq$RPKM[-(1:29)])
colnames(cancerRpkmValues) <- substring(cancerFileNames[1], 14, 27)
rownames(cancerRpkmValues) <- genes
for (i in 2:nCancer){
  cat("Reading cancer rpkm values from file", i, "of", nCancer, "\n")
  tmpRnaseq <- read.table(paste(dataFileDir, cancerFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  cancerRpkmValues[[substring(cancerFileNames[i], 14, 27)]] <- tmpRnaseq$RPKM[-(1:29)]
}

# Transpose and filter out bad data
ggNormal <- goodGenes(t(normalRpkmValues), verbose=3)
ggCancer <- goodGenes(t(cancerRpkmValues), verbose=3)
normalRnaseq <- as.data.frame(t(normalRpkmValues)[, ggNormal & ggCancer])
cancerRnaseq <- as.data.frame(t(cancerRpkmValues)[, ggNormal & ggCancer])
genesList <- colnames(normalRnaseq)

save(normalRnaseq, file="../Rdata/normalRnaseqAllGenes.Rdata")
save(cancerRnaseq, file="../Rdata/cancerRnaseqAllGenes.Rdata")
save(genesList,    file="../Rdata/genesList.Rdata")
