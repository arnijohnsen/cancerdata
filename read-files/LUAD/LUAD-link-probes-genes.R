# Read annotation file and load probe and gene list
cat("Read probe annotation, probe list and gene list\n")
probe.annotation <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", 
                               header=TRUE, sep="\t", quote="\"", stringsAsFactors=FALSE)
load("../Rdata/LUAD/info/LUAD-probes.Rdata")
load("../Rdata/LUAD/info/LUAD-genes.Rdata")

# Filter probe.annotation with probes in use
probe.annotation <- probe.annotation[probe.annotation$TargetID %in% LUAD.probes,]
probe.gene.info <- data.frame(probe = probe.annotation$TargetID, 
                            genes = I(strsplit(probe.annotation$UCSC_REFGENE_NAME, ";")), 
			    types = I(strsplit(probe.annotation$UCSC_REFGENE_GROUP, ";")))

# Create vectors to link probes to genes
cat("Running loop to link probes to genes\n")
probe.loop <- character(0)
genes.loop <- character(0)

# Run loop to generate a linked probe-gene list
# Every loop iteration adds one probe and all associated genes to the list 
# (which can result in multiple entries for one probe if it's associated with
#  more than one gene)
l <- 0
n <- dim(probe.gene.info)[1]
pb <- txtProgressBar(min=1, max=n, style=3)
for (i in 1:n){
  setTxtProgressBar(pb, i)
  probe <- as.character(probe.gene.info[i,1])
  genes <- as.vector(probe.gene.info[i,2][[1]])
  types <- as.vector(probe.gene.info[i,3][[1]])
  # Check if only one gene
  if (length(unique(genes)) == 1){
    probe.loop <- c(probe.loop, probe)
    genes.loop <- c(genes.loop, unique(genes))
    l <- l+1
  } else {
    # Multiple genes, select only promoter associated genes
    genes <- genes[grep("TSS200|5'UTR", types)]
    genes <- unique(genes)
    for (k in 1:length(genes)){
      probe.loop <- c(probe.loop, probe)
      genes.loop <- c(genes.loop, genes[k])
      l <- l+1
    }
  }
  if (length(genes.loop) != length(probe.loop)){
    cat("error at i =", i, "\n")
    break
  }
}
cat("\n")

# After list is generated, remove all entries with genes which aren't in 
# the list of genes use 
idx <- genes.loop %in% LUAD.genes

# Save data to file and exit
LUAD.linked.probes.genes <- data.frame(probes = probe.loop[idx], genes = genes.loop[idx])
save(LUAD.linked.probes.genes, file="../Rdata/LUAD/info/LUAD-linked-probes-genes.Rdata")
quit(save="no")
