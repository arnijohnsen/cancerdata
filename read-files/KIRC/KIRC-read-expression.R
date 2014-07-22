# Load libraries
library(WGCNA)

# Set data file directory and get information about file names and barcodes
cat("Retrieving file list\n")
data.file.dir = "../rawdata/KIRC/RNASeq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
file.sample.map <- read.table("../rawdata/KIRC/RNASeq/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
colnames(file.sample.map) <- c("filename", "barcode")
file.sample.map <- file.sample.map[grep("genes.normalized_results", file.sample.map$filename),]
# Use only first 14 letters of barcode
file.sample.map$barcode <- substring(gsub(".*,TCGA", "TCGA", file.sample.map$barcode), 1, 14)

n <- length(file.sample.map$filename)

# Read first file to get gene names and dimension of data frame
tmp.rnaseq <- read.table(paste(data.file.dir, file.sample.map$filename[1], sep=""), 
                        header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
# Create list of genes and fix one wrong gene name
genes <- gsub("\\|.*", "", tmp.rnaseq$gene_id[-(1:29)])
genes[grep("SLC35E2", genes)][2] <- "SLC35E2B" # Fix one wrong gene name

# Create data frames for beta values
normal.rpkm.values <- data.frame(rownames = genes)
cancer.rpkm.values <- data.frame(rownames = genes)

cat("Reading cancer rpkm values\n")
pb <- txtProgressBar(min=1, max=n, style=3)
contains.normal <- FALSE
# In each loop iteration, one file is read and from its barcode it's determined
#  wether it's normal or cancer, and that data is saved to the corresponding frame
for (i in 1:n){
  setTxtProgressBar(pb, i)
  tmp.rnaseq <- read.table(paste(data.file.dir, file.sample.map$filename[i], sep=""), 
                          header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
  if(substring(file.sample.map$barcode[i], 14,14) == "0"){
    cancer.rpkm.values[[file.sample.map$barcode[i]]] <- tmp.rnaseq$normalized_count[-(1:29)]
  }else{
    normal.rpkm.values[[file.sample.map$barcode[i]]] <- tmp.rnaseq$normalized_count[-(1:29)]
    contains.normal <- TRUE
  }
}
cat("\n")
# Fix rownames
cancer.rpkm.values$rownames <- NULL
rownames(cancer.rpkm.values) <- genes
if(contains.normal){
  normal.rpkm.values$rownames <- NULL
  rownames(normal.rpkm.values) <- genes
}

# Transpose and filter out bad data
gg.cancer <- goodGenes(t(cancer.rpkm.values), verbose=3)
if(contains.normal){
  gg.normal <- goodGenes(t(normal.rpkm.values), verbose=3)
}
# Save files, two versions depending on if there are any normal samples
if(contains.normal){
  KIRC.CEA <- as.data.frame(t(cancer.rpkm.values)[,gg.normal & gg.cancer])
  KIRC.NEA <- as.data.frame(t(normal.rpkm.values)[,gg.normal & gg.cancer])
  KIRC.genes <- colnames(KIRC.CEA)
  KIRC.CEA.samples <- rownames(KIRC.CEA)
  KIRC.NEA.samples <- rownames(KIRC.NEA)
  save(KIRC.CEA,         file="../Rdata/KIRC/data/KIRC-CEA.Rdata")
  save(KIRC.NEA,         file="../Rdata/KIRC/data/KIRC-NEA.Rdata")
  save(KIRC.genes,      file="../Rdata/KIRC/info/KIRC-genes.Rdata")
  save(KIRC.CEA.samples, file="../Rdata/KIRC/info/KIRC-CEA-samples.Rdata")
  save(KIRC.NEA.samples, file="../Rdata/KIRC/info/KIRC-NEA-samples.Rdata")
}else{
  KIRC.CEA <- as.data.frame(t(cancer.rpkm.values)[,gg.cancer])
  KIRC.genes <- colnames(KIRC.CEA)
  KIRC.CEA.samples <- rownames(KIRC.CEA)
  save(KIRC.CEA,         file="../Rdata/KIRC/data/KIRC-CEA.Rdata")
  save(KIRC.genes,      file="../Rdata/KIRC/info/KIRC-genes.Rdata")
  save(KIRC.CEA.samples, file="../Rdata/KIRC/info/KIRC-CEA-samples.Rdata")
}
quit(save="no")
