# Load libraries
library(WGCNA)

# Read normal and cancer tissue rna sequencing files (regex matches only normal/cancer files)
cat("Retrieving file list\n")
data.file.dir = "../rawdata/KIRC/RNASeq/RNASeq/UNC__IlluminaHiSeq_RNASeq/Level_3/"
normal.file.names <- list.files(data.file.dir, pattern=".*1[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
cancer.file.names <- list.files(data.file.dir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
# Rna sequencing files contain multiple files for each barcode, select only gene expression
normal.file.names <- normal.file.names[grep("gene", normal.file.names)]
cancer.file.names <- cancer.file.names[grep("gene", cancer.file.names)]

n.normal <- length(normal.file.names)
# Read first file to get gene names and dimension of data frame
cat("Reading normal rpkm values from file", 1, "of", n.normal, "\n")
tmp.rnaseq <- read.table(paste(data.file.dir, normal.file.names[1], sep=""), 
                        header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmp.rnaseq$gene[-(1:29)])
genes[16272] <- "SLC35E2B" # Fix one wrong gene name

# Create data frame, don't use first 29 genes as they don't have names
normal.rpkm.values <- data.frame(x = tmp.rnaseq$RPKM[-(1:29)])
colnames(normal.rpkm.values) <- substring(gsub(".*TCGA", "TCGA", normal.file.names[1]), 1, 14)
rownames(normal.rpkm.values) <- genes

# Run loop over the rest of normal files
cat("Reading normal rpkm values\n")
pb <- txtProgressBar(min=1, max=n.normal, style=3)
for (i in 2:n.normal){
  setTxtProgressBar(pb, i)
  tmp.rnaseq <- read.table(paste(data.file.dir, normal.file.names[i], sep=""), 
                          header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  normal.rpkm.values[[substring(gsub(".*TCGA", "TCGA", normal.file.names[i]), 1, 14)]] <- tmp.rnaseq$RPKM[-(1:29)]
}
cat("\n")

n.cancer <- length(cancer.file.names)
# Read first file to get gene names and dimension of data frame
cat("Reading cancer rpkm values from file", 1, "of", n.cancer, "\n")
tmp.rnaseq <- read.table(paste(data.file.dir, cancer.file.names[1], sep=""), 
                        header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmp.rnaseq$gene[-(1:29)])
genes[16272] <- "SLC35E2B" # Fix one wrong gene name

# Create data frame, don't use first 29 genes as they don't have names
cancer.rpkm.values <- data.frame(x = tmp.rnaseq$RPKM[-(1:29)])
colnames(cancer.rpkm.values) <- substring(gsub(".*TCGA", "TCGA", cancer.file.names[1]), 1, 14)
rownames(cancer.rpkm.values) <- genes

# Run loop over the rest of normal files
cat("Reading cancer rpkm values\n")
pb <- txtProgressBar(min=1, max=n.cancer, style=3)
for (i in 2:n.cancer){
  setTxtProgressBar(pb, i)
  tmp.rnaseq <- read.table(paste(data.file.dir, cancer.file.names[i], sep=""), 
                          header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  cancer.rpkm.values[[substring(gsub(".*TCGA", "TCGA", cancer.file.names[i]), 1, 14)]] <- tmp.rnaseq$RPKM[-(1:29)]
}
cat("\n")

# Transpose and filter out bad data
gg.normal <- goodGenes(t(normal.rpkm.values), verbose=3)
gg.cancer <- goodGenes(t(cancer.rpkm.values), verbose=3)
KIRC.NEA <- as.data.frame(t(normal.rpkm.values)[, gg.normal & gg.cancer])
KIRC.CEA <- as.data.frame(t(cancer.rpkm.values)[, gg.normal & gg.cancer])
KIRC.genes <- colnames(KIRC.NEA)
KIRC.NEA.samples <- rownames(KIRC.NEA)
KIRC.CEA.samples <- rownames(KIRC.CEA)

# Save data to file and exit
cat("Saving data to file\n")
save(KIRC.NEA,         file="../Rdata/KIRC/data/KIRC-NEA.Rdata")
save(KIRC.CEA,         file="../Rdata/KIRC/data/KIRC-CEA.Rdata")
save(KIRC.genes,       file="../Rdata/KIRC/info/KIRC-genes.Rdata")
save(KIRC.NEA.samples, file="../Rdata/KIRC/info/KIRC-NEA-samples.Rdata")
save(KIRC.CEA.samples, file="../Rdata/KIRC/info/KIRC-CEA-samples.Rdata")
quit(save="no")
