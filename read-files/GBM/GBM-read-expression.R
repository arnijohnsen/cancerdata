# Load libraries
library(WGCNA)

# Read cancer tissue rna sequencing files (regex matches only cancer files)
cat("Retrieving file list\n")
data.file.dir = "../rawdata/GBM/RNASeq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
cancer.file.names <- list.files(data.file.dir)
# Get only annotated gene files
cancer.file.names <- cancer.file.names[grep("genes.normalized_results", cancer.file.names)]

# File names don't include barcodes, read FILE_SAMPLE_MAP to get barcodes
link <- read.table("../rawdata/GBM/RNASeq/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
# Use only first 14 letters of barcode
link$barcode.s. <- substring(gsub(".*,TCGA", "TCGA", link$barcode.s.), 1, 14)

n.cancer <- length(cancer.file.names)
# Read first file to get gene names and dimension of data frame
cat("Reading cancer rpkm values from file", 1, "of", n.cancer, "\n")
tmp.rnaseq <- read.table(paste(data.file.dir, cancer.file.names[1], sep=""), 
                        header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmp.rnaseq$gene_id[-(1:29)])
genes[grep("SLC35E2", genes)][2] <- "SLC35E2B" # Fix one wrong gene name

# Create data frame, don't use first 29 genes as they don't have names
cancer.rpkm.values <- data.frame(x = tmp.rnaseq$normalized_count[-(1:29)])
colnames(cancer.rpkm.values) <- link$barcode.s.[link$filename == cancer.file.names[1]]
rownames(cancer.rpkm.values) <- genes
cat("Reading cancer rpkm values\n")
pb <- txtProgressBar(min=1, max=n.cancer, style=3)
for (i in 2:n.cancer){
  setTxtProgressBar(pb, i)
  tmp.rnaseq <- read.table(paste(data.file.dir, cancer.file.names[i], sep=""), 
                          header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  cancer.rpkm.values[[ link$barcode.s.[link$filename == cancer.file.names[i]] ]] <- tmp.rnaseq$normalized_count[-(1:29)]
}
cat("\n")
# Transpose and filter out bad data
gg.cancer <- goodGenes(t(cancer.rpkm.values), verbose=3)
GBM.CEA <- as.data.frame(t(cancer.rpkm.values)[,gg.cancer])
GBM.genes <- colnames(GBM.CEA)
GBM.CEA.samples <- rownames(GBM.CEA)

# Save data to file and exit
cat("Saving data to file\n")
save(GBM.CEA,         file="../Rdata/GBM/data/GBM-CEA.Rdata")
save(GBM.genes,       file="../Rdata/GBM/info/GBM-genes.Rdata")
save(GBM.CEA.samples, file="../Rdata/GBM/info/GBM-CEA-samples.Rdata")
quit(save="no")
