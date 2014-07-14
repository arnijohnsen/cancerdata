# Load libraries
library(WGCNA)

# Get info about promoter probes, defined as TSS200 or 5'UTR probes
cat("Reading probes annotation file\n")
probe.annotation <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"")
idx <- grep("TSS200|5'UTR", probe.annotation$UCSC_REFGENE_GROUP)
prom.probes <- probe.annotation$TargetID[idx]

# Read cancer tissue methylation files (regex matches only cancer files)
cat("Retrieving file list\n")
data.file.dir <- "../rawdata/KIRC/Methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/"
normal.file.names <- list.files(data.file.dir, pattern=".*1[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")
cancer.file.names <- list.files(data.file.dir, pattern=".*0[A-Z0-9]{2}-[A-Z0-9]{3}-[A-Z0-9]{4}-[0-9]{2}.*")

# Create new data.frame containing beta values, 
# with each normal sample in one column and each probe in one row
normal.beta.values <- data.frame(rownames = prom.probes)
n.normal <- length(normal.file.names)
cat("Reading normal beta values\n")
pb <- txtProgressBar(min=1, max=n.normal, style=3)
for (i in 1:n.normal){
  setTxtProgressBar(pb, i)
  tmp.methyl <- read.table(paste(data.file.dir, normal.file.names[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmp.methyl[tmp.methyl$Composite.Element.REF %in% prom.probes, ]
  normal.beta.values[[substring(gsub(".*TCGA", "TCGA", normal.file.names[i]),1,14)]] <- tmp$Beta_value 
}
cat("\n")
normal.beta.values$rownames <- NULL
rownames(normal.beta.values) <- prom.probes

# Create new data.frame containing beta values, 
# with each cancer sample in one column and each probe in one row
cancer.beta.values <- data.frame(rownames = prom.probes)
n.cancer <- length(cancer.file.names)
cat("Reading cancer beta values\n")
pb <- txtProgressBar(min=1, max=n.cancer, style=3)
for (i in 1:n.cancer){
  setTxtProgressBar(pb, i)
  tmp.methyl <- read.table(paste(data.file.dir, cancer.file.names[i], sep=""), header=TRUE, sep="\t", quote="\"", skip=1)
  tmp <- tmp.methyl[tmp.methyl$Composite.Element.REF %in% prom.probes, ]
  cancer.beta.values[[substring(gsub(".*TCGA", "TCGA", cancer.file.names[i]),1,14)]] <- tmp$Beta_value 
}
cat("\n")
cancer.beta.values$rownames <- NULL
rownames(cancer.beta.values) <- prom.probes

# Transpose and filter out bad data
gg.cancer <- goodGenes(t(cancer.beta.values), verbose=3)

KIRC.CMP <- as.data.frame(t(cancer.beta.values)[,gg.cancer])
KIRC.probes <- colnames(KIRC.CMP)
KIRC.CMP.samples <- rownames(KIRC.CMP)

# Save to files
save(KIRC.CMP,         file="../Rdata/KIRC/data/KIRC-CMP.Rdata")
save(KIRC.probes,      file="../Rdata/KIRC/info/KIRC-probes.Rdata")
save(KIRC.CMP.samples, file="../Rdata/KIRC/info/KIRC-CMP-samples.Rdata")
quit(save="no")
