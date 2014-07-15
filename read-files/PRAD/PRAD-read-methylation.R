# Load libraries
library(WGCNA)

# Get info about promoter probes, defined as TSS200 or 5'UTR probes
cat("Reading probes annotation file\n")
probe.annotation <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"")
idx <- grep("TSS200|5'UTR", probe.annotation$UCSC_REFGENE_GROUP)
prom.probes <- probe.annotation$TargetID[idx]

# Read cancer tissue methylation files (regex matches only cancer files)
cat("Retrieving file list\n")
data.file.dir = "../rawdata/PRAD/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
file.sample.map <- read.table("../rawdata/PRAD/Methylation/FILE_SAMPLE_MAP.txt", header=T, sep="\t")
colnames(file.sample.map) <- c("filename", "barcode")
# Use only first 14 letters of barcode
file.sample.map$barcode <- substring(gsub(".*,TCGA", "TCGA", file.sample.map$barcode), 1, 14)

normal.beta.values <- data.frame(rownmaes = prom.probes)
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
  PRAD.CMP <- as.data.frame(t(cancer.beta.values)[,gg.normal & gg.cancer])
  PRAD.NMP <- as.data.frame(t(normal.beta.values)[,gg.normal & gg.cancer])
  PRAD.probes <- colnames(PRAD.CMP)
  PRAD.CMP.samples <- rownames(PRAD.CMP)
  PRAD.NMP.samples <- rownames(PRAD.NMP)
  save(PRAD.CMP,         file="../Rdata/PRAD/data/PRAD-CMP.Rdata")
  save(PRAD.NMP,         file="../Rdata/PRAD/data/PRAD-NMP.Rdata")
  save(PRAD.probes,      file="../Rdata/PRAD/info/PRAD-probes.Rdata")
  save(PRAD.CMP.samples, file="../Rdata/PRAD/info/PRAD-CMP-samples.Rdata")
  save(PRAD.NMP.samples, file="../Rdata/PRAD/info/PRAD-NMP-samples.Rdata")
}else{
  PRAD.CMP <- as.data.frame(t(cancer.beta.values)[,gg.cancer])
  PRAD.probes <- colnames(PRAD.CMP)
  PRAD.CMP.samples <- rownames(PRAD.CMP)
  save(PRAD.CMP,         file="../Rdata/PRAD/data/PRAD-CMP.Rdata")
  save(PRAD.probes,      file="../Rdata/PRAD/info/PRAD-probes.Rdata")
  save(PRAD.CMP.samples, file="../Rdata/PRAD/info/PRAD-CMP-samples.Rdata")
}
quit(save="no")
