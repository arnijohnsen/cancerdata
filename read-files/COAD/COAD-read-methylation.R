# Load libraries
library(WGCNA)

# Get info about promoter probes, defined as TSS200 or 5'UTR probes
cat("Reading probes annotation file\n")
probe.annotation <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"")
idx <- grep("TSS200|5'UTR", probe.annotation$UCSC_REFGENE_GROUP)
prom.probes <- probe.annotation$TargetID[idx]

# Set data file directory and get information about file names and barcodes
cat("Retrieving file list\n")
data.file.dir = "../rawdata/COAD/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
file.sample.map <- read.table("../rawdata/COAD/Methylation/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
colnames(file.sample.map) <- c("filename", "barcode")
# Use only first 14 letters of barcode
file.sample.map$barcode <- substring(gsub(".*,TCGA", "TCGA", file.sample.map$barcode), 1, 14)

# Create data frames for beta values
normal.beta.values <- data.frame(rownames = prom.probes)
cancer.beta.values <- data.frame(rownames = prom.probes)
n <- length(file.sample.map$filename)
cat("Reading beta values\n")
pb <- txtProgressBar(min=1, max=n, style=3)
contains.normal <- FALSE
# In each loop iteration, one file is read and from its barcode it's determined
#  wether it's normal or cancer, and that data is saved to the corresponding frame
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
# Fix rownames
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
# Save files, two versions depending on if there are any normal samples
if(contains.normal){
  COAD.CMP <- as.data.frame(t(cancer.beta.values)[,gg.normal & gg.cancer])
  COAD.NMP <- as.data.frame(t(normal.beta.values)[,gg.normal & gg.cancer])
  COAD.probes <- colnames(COAD.CMP)
  COAD.CMP.samples <- rownames(COAD.CMP)
  COAD.NMP.samples <- rownames(COAD.NMP)
  save(COAD.CMP,         file="../Rdata/COAD/data/COAD-CMP.Rdata")
  save(COAD.NMP,         file="../Rdata/COAD/data/COAD-NMP.Rdata")
  save(COAD.probes,      file="../Rdata/COAD/info/COAD-probes.Rdata")
  save(COAD.CMP.samples, file="../Rdata/COAD/info/COAD-CMP-samples.Rdata")
  save(COAD.NMP.samples, file="../Rdata/COAD/info/COAD-NMP-samples.Rdata")
}else{
  COAD.CMP <- as.data.frame(t(cancer.beta.values)[,gg.cancer])
  COAD.probes <- colnames(COAD.CMP)
  COAD.CMP.samples <- rownames(COAD.CMP)
  save(COAD.CMP,         file="../Rdata/COAD/data/COAD-CMP.Rdata")
  save(COAD.probes,      file="../Rdata/COAD/info/COAD-probes.Rdata")
  save(COAD.CMP.samples, file="../Rdata/COAD/info/COAD-CMP-samples.Rdata")
}
quit(save="no")
