# Load libraries
library(WGCNA)

# Read normal and cancer tissue rna sequencing files (regex matches only normal/cancer files)
cat("Retrieving file list\n")
dataFileDir = "../rawdata/KIRC/RNASeq/RNASeq/UNC__IlluminaHiSeq_RNASeq/Level_3/"
normalFileNames <- list.files(dataFileDir, pattern=".*1[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
cancerFileNames <- list.files(dataFileDir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
# Rna sequencing files contain multiple files for each barcode, select only gene expression
normalFileNames <- normalFileNames[grep("gene", normalFileNames)]
cancerFileNames <- cancerFileNames[grep("gene", cancerFileNames)]

nNormal <- length(normalFileNames)
# Read first file to get gene names and dimension of data frame
cat("Reading normal rpkm values from file", 1, "of", nNormal, "\n")
tmpRnaseq <- read.table(paste(dataFileDir, normalFileNames[1], sep=""), 
                        header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmpRnaseq$gene[-(1:29)])
genes[16272] <- "SLC35E2B" # Fix one wrong gene name

# Create data frame, don't use first 29 genes as they don't have names
normalRpkmValues <- data.frame(x = tmpRnaseq$RPKM[-(1:29)])
colnames(normalRpkmValues) <- substring(gsub(".*TCGA", "TCGA", normalFileNames[1]), 1, 14)
rownames(normalRpkmValues) <- genes

# Run loop over the rest of normal files
cat("Reading normal rpkm values\n")
pb <- txtProgressBar(min=1, max=nNormal, style=3)
for (i in 2:nNormal){
  setTxtProgressBar(pb, i)
  tmpRnaseq <- read.table(paste(dataFileDir, normalFileNames[i], sep=""), 
                          header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  normalRpkmValues[[substring(gsub(".*TCGA", "TCGA", normalFileNames[i]), 1, 14)]] <- tmpRnaseq$RPKM[-(1:29)]
}
cat("\n")

nCancer <- length(cancerFileNames)
# Read first file to get gene names and dimension of data frame
cat("Reading cancer rpkm values from file", 1, "of", nCancer, "\n")
tmpRnaseq <- read.table(paste(dataFileDir, cancerFileNames[1], sep=""), 
                        header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmpRnaseq$gene[-(1:29)])
genes[16272] <- "SLC35E2B" # Fix one wrong gene name

# Create data frame, don't use first 29 genes as they don't have names
cancerRpkmValues <- data.frame(x = tmpRnaseq$RPKM[-(1:29)])
colnames(cancerRpkmValues) <- substring(gsub(".*TCGA", "TCGA", cancerFileNames[1]), 1, 14)
rownames(cancerRpkmValues) <- genes

# Run loop over the rest of normal files
cat("Reading cancer rpkm values\n")
pb <- txtProgressBar(min=1, max=nNormal, style=3)
for (i in 2:nCancer){
  setTxtProgressBar(pb, i)
  tmpRnaseq <- read.table(paste(dataFileDir, cancerFileNames[i], sep=""), 
                          header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  cancerRpkmValues[[substring(gsub(".*TCGA", "TCGA", cancerFileNames[i]), 1, 14)]] <- tmpRnaseq$RPKM[-(1:29)]
}
cat("\n")

# Transpose and filter out bad data
ggNormal <- goodGenes(t(normalRpkmValues), verbose=3)
ggCancer <- goodGenes(t(cancerRpkmValues), verbose=3)
normalRnaseq <- as.data.frame(t(normalRpkmValues)[, ggNormal & ggCancer])
cancerRnaseq <- as.data.frame(t(cancerRpkmValues)[, ggNormal & ggCancer])
genesList <- colnames(normalRnaseq)
normalRnaseqSampList <- rownames(normalRnaseq)
cancerRnaseqSampList <- rownames(cancerRnaseq)

# Save data to file and exit
cat("Saving data to file\n")
save(normalRnaseq, file="../Rdata/KIRC/data/KIRC-NEA.Rdata")
save(cancerRnaseq, file="../Rdata/KIRC/data/KIRC-CEA.Rdata")
save(genesList,    file="../Rdata/KIRC/info/genesList.Rdata")
save(normalRnaseqSampList, cancerRnaseqSampList, file="../Rdata/KIRC/info/rnaseqSampList.Rdata")
quit(save="no")
