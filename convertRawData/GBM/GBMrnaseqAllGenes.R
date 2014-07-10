# Load libraries
library(WGCNA)

# Read cancer tissue rna sequencing files (regex matches only cancer files)
cat("Retrieving file list\n")
dataFileDir = "../rawdata/GBM/RNASeq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
cancerFileNames <- list.files(dataFileDir)
# Get only annotated gene files
cancerFileNames <- cancerFileNames[grep("genes.normalized_results", cancerFileNames)]

link <- read.table("../rawdata/GBM/RNASeq/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
link$barcode.s. <- substring(gsub(".*,TCGA", "TCGA", link$barcode.s.), 1, 14)

nCancer <- length(cancerFileNames)
# Read first file to get gene names and dimension of data frame
cat("Reading cancer rpkm values from file", 1, "of", nCancer, "\n")
tmpRnaseq <- read.table(paste(dataFileDir, cancerFileNames[1], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmpRnaseq$gene_id[-(1:29)])
# Fix one wrong names
genes[grep("SLC35E2", genes)][2] <- "SLC35E2B"
cancerRpkmValues <- data.frame(x = tmpRnaseq$normalized_count[-(1:29)])
colnames(cancerRpkmValues) <- link$barcode.s.[link$filename == cancerFileNames[1]]
rownames(cancerRpkmValues) <- genes
cat("Reading cancer rpkm values\n")
pb <- txtProgressBar(min=1, max=nCancer, style=3)
for (i in 2:nCancer){
  setTxtProgressBar(pb, i)
  tmpRnaseq <- read.table(paste(dataFileDir, cancerFileNames[i], sep=""), header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  cancerRpkmValues[[ link$barcode.s.[link$filename == cancerFileNames[i]] ]] <- tmpRnaseq$normalized_count[-(1:29)]
}
cat("\n")
print(dim(cancerRpkmValues))
print(nCancer)
print(colnames(cancerRpkmValues))
# Transpose and filter out bad data
ggCancer <- goodGenes(t(cancerRpkmValues), verbose=3)
cancerRnaseq <- as.data.frame(t(cancerRpkmValues)[,ggCancer])
genesList <- colnames(cancerRnaseq)
cancerRnaseqSampList <- rownames(cancerRnaseq)

cat("Saving data to file\n")
save(cancerRnaseq, file="../Rdata/GBM/data/GBM-CEA.Rdata")
save(genesList,    file="../Rdata/GBM/info/genesList.Rdata")
save(cancerRnaseqSampList, file="../Rdata/GBM/info/rnaseqSampList.Rdata")
quit(save="no")
