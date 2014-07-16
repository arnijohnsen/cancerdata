# Load libraries
library(WGCNA)

# Get info about promoter probes, defined as TSS200 or 5'UTR probes
cat("Reading probes annotation file\n")
probe.annotation <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"")
idx <- grep("TSS200|5'UTR", probe.annotation$UCSC_REFGENE_GROUP)
prom.probes <- probe.annotation$TargetID[idx]

# Read cancer tissue methylation files (regex matches only cancer files)
cat("Retrieving file list\n")
data.file.dir = "../rawdata/LIHC/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
file.sample.map <- read.table("../rawdata/LIHC/Methylation/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
colnames(file.sample.map) <- c("filename", "barcode")
# Use only first 14 letters of barcode
file.sample.map$barcode <- substring(gsub(".*,TCGA", "TCGA", file.sample.map$barcode), 1, 14)

normal.beta.values <- data.frame(rownames = prom.probes)
cancer.beta.values <- data.frame(rownames = prom.probes)
n <- length(file.sample.map$filename)
cat("Reading beta values\n")
pb <- txtProgressBar(min=1, max=n, style=3)
contains.normal <- FALSE
for (i in 1:n){
  setTxtProgressBar(pb, i)
  tmp.methyl <- read.table(paste(data.file.dir, file.sample.map$filename[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmp.methyl[tmp.methyl$Composite.Element.REF %in% prom.probes, ]
  if(substring(file.sample.map$barcode[i], 14,14) == "0"){
    cancer.beta.values[[file.sample.map$barcode[i]]] <- tmp$Beta_value 
  }else{
    normal.beta.values[[file.sample.map$barcode[i]]] <- tmp$Beta_value 
    contains.normal <- TRUE
  }
}
cat("\n")
cancer.beta.values$rownames <- NULL
rownames(cancer.beta.values) <- prom.probes
if(contains.normal){
  normal.beta.values$rownames <- NULL
  rownames(normal.beta.values) <- prom.probes
}

# Transpose and filter out bad data
gg.cancer <- goodGenes(t(cancer.beta.values), verbose=3)
if(contains.normal){
  gg.normal <- goodGenes(t(normal.beta.values), verbose=3)
}
if(contains.normal){
  LIHC.CMP <- as.data.frame(t(cancer.beta.values)[,gg.normal & gg.cancer])
  LIHC.NMP <- as.data.frame(t(normal.beta.values)[,gg.normal & gg.cancer])
  LIHC.probes <- colnames(LIHC.CMP)
  LIHC.CMP.samples <- rownames(LIHC.CMP)
  LIHC.NMP.samples <- rownames(LIHC.NMP)
  save(LIHC.CMP,         file="../Rdata/LIHC/data/LIHC-CMP.Rdata")
  save(LIHC.NMP,         file="../Rdata/LIHC/data/LIHC-NMP.Rdata")
  save(LIHC.probes,      file="../Rdata/LIHC/info/LIHC-probes.Rdata")
  save(LIHC.CMP.samples, file="../Rdata/LIHC/info/LIHC-CMP-samples.Rdata")
  save(LIHC.NMP.samples, file="../Rdata/LIHC/info/LIHC-NMP-samples.Rdata")
}else{
  LIHC.CMP <- as.data.frame(t(cancer.beta.values)[,gg.cancer])
  LIHC.probes <- colnames(LIHC.CMP)
  LIHC.CMP.samples <- rownames(LIHC.CMP)
  save(LIHC.CMP,         file="../Rdata/LIHC/data/LIHC-CMP.Rdata")
  save(LIHC.probes,      file="../Rdata/LIHC/info/LIHC-probes.Rdata")
  save(LIHC.CMP.samples, file="../Rdata/LIHC/info/LIHC-CMP-samples.Rdata")
}
quit(save="no")
