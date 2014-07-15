# Load libraries
library(WGCNA)

# Read cancer tissue methylation files (regex matches only cancer files)
cat("Retrieving file list\n")
data.file.dir = "../rawdata/LIHC/RNASeq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
file.sample.map <- read.table("../rawdata/LIHC/RNASeq/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
colnames(file.sample.map) <- c("filename", "barcode")
file.sample.map <- file.sample.map[grep("genes.normalized_results", file.sample.map$filename),]
# Use only first 14 letters of barcode
file.sample.map$barcode <- substring(gsub(".*,TCGA", "TCGA", file.sample.map$barcode), 1, 14)

n <- length(file.sample.map$filename)

# Read first file to get gene names and dimension of data frame
tmp.rnaseq <- read.table(paste(data.file.dir, file.sample.map$filename[1], sep=""), 
                        header=TRUE, sep="\t", quote="\"", stringsAsFactors = FALSE)
genes <- gsub("\\|.*", "", tmp.rnaseq$gene_id[-(1:29)])
genes[grep("SLC35E2", genes)][2] <- "SLC35E2B" # Fix one wrong gene name
normal.rpkm.values <- data.frame(rownames = genes)
cancer.rpkm.values <- data.frame(rownames = genes)

cat("Reading cancer rpkm values\n")
pb <- txtProgressBar(min=1, max=n, style=3)
contains.normal <- FALSE
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
if(contains.normal){
  LIHC.CEA <- as.data.frame(t(cancer.rpkm.values)[,gg.normal & gg.cancer])
  LIHC.NEA <- as.data.frame(t(normal.rpkm.values)[,gg.normal & gg.cancer])
  LIHC.genes <- colnames(LIHC.CEA)
  LIHC.CEA.samples <- rownames(LIHC.CEA)
  LIHC.NEA.samples <- rownames(LIHC.NEA)
  save(LIHC.CEA,         file="../Rdata/LIHC/data/LIHC-CEA.Rdata")
  save(LIHC.NEA,         file="../Rdata/LIHC/data/LIHC-NEA.Rdata")
  save(LIHC.genes,      file="../Rdata/LIHC/info/LIHC-genes.Rdata")
  save(LIHC.CEA.samples, file="../Rdata/LIHC/info/LIHC-CEA-samples.Rdata")
  save(LIHC.NEA.samples, file="../Rdata/LIHC/info/LIHC-NEA-samples.Rdata")
}else{
  LIHC.CEA <- as.data.frame(t(cancer.rpkm.values)[,gg.cancer])
  LIHC.genes <- colnames(LIHC.CEA)
  LIHC.CEA.samples <- rownames(LIHC.CEA)
  save(LIHC.CEA,         file="../Rdata/LIHC/data/LIHC-CEA.Rdata")
  save(LIHC.genes,      file="../Rdata/LIHC/info/LIHC-genes.Rdata")
  save(LIHC.CEA.samples, file="../Rdata/LIHC/info/LIHC-CEA-samples.Rdata")
}
quit(save="no")
