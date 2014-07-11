# Read annotation file and load probe and gene list
cat("Read probe annotation, probe list and gene list\n")
probeAnnotations <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", 
                               header=TRUE, sep="\t", quote="\"", stringsAsFactors=FALSE)
load("../Rdata/BRCA/info/probeList.Rdata")
load("../Rdata/BRCA/info/genesList.Rdata")

# Filter probeAnnotations with probes in use
probeAnnotations <- probeAnnotations[probeAnnotations$TargetID %in% probeList,]
probeGeneInfo <- data.frame(probe = probeAnnotations$TargetID, 
                            genes = I(strsplit(probeAnnotations$UCSC_REFGENE_NAME, ";")), 
			    types = I(strsplit(probeAnnotations$UCSC_REFGENE_GROUP, ";")))

# Create vectors to link probes to genes
cat("Running loop to link probes to genes\n")
probeLoop <- character(0)
genesLoop <- character(0)

# Run loop to generate a linked probe-gene list
# Every loop iteration adds one probe and all associated genes to the list 
# (which can result in multiple entries for one probe if it's associated with
#  more than one gene)
l <- 0
n <- dim(probeGeneInfo)[1]
pb <- txtProgressBar(min=1, max=n, style=3)
for (i in 1:n){
  setTxtProgressBar(pb, i)
  probe <- as.character(probeGeneInfo[i,1])
  genes <- as.vector(probeGeneInfo[i,2][[1]])
  types <- as.vector(probeGeneInfo[i,3][[1]])
  # Check if only one gene
  if (length(unique(genes)) == 1){
    probeLoop <- c(probeLoop, probe)
    genesLoop <- c(genesLoop, unique(genes))
    l <- l+1
  } else {
    # Multiple genes, select only promoter associated genes
    genes <- genes[grep("TSS200|5'UTR", types)]
    genes <- unique(genes)
    for (k in 1:length(genes)){
      probeLoop <- c(probeLoop, probe)
      genesLoop <- c(genesLoop, genes[k])
      l <- l+1
    }
  }
  if (length(genesLoop) != length(probeLoop)){
    cat("error at i =", i, "\n")
    break
  }
}
cat("\n")

# After list is generated, remove all entries with genes which aren't in 
# the list of genes use 
idx <- genesLoop %in% genesList

# Save data to file and exit
linkedProbesGenes <- data.frame(probes = probeLoop[idx], genes = genesLoop[idx])
save(linkedProbesGenes, file="../Rdata/BRCA/info/linkedProbesGenes.Rdata")
quit(save="no")
