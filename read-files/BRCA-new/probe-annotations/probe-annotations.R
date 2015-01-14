cat("Reading probes annotation file..\n")
prb.ant.raw <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"", stringsAsFactors=FALSE)

cat("Setting probe status..\n")
prb.ant.raw$status <- 0

cat("  ..gene bodies\n")
is.body   <- !grepl("TSS200|TSS1500|1stExon|5'UTR|3'UTR|^$", prb.ant.raw$UCSC_REFGENE_GROUP)
is.island <-  grepl("Island", prb.ant.raw$RELATION_TO_UCSC_CPG_ISLAND)
is.shore  <-  grepl("Shore",  prb.ant.raw$RELATION_TO_UCSC_CPG_ISLAND)
is.none   <- !is.island & !is.shore
prb.ant.raw$status[is.body & is.island] <- "body.island"
prb.ant.raw$status[is.body & is.shore]  <- "body.shore"
prb.ant.raw$status[is.body & is.none]   <- "body.none"

cat("  ..enhancers\n")
prb.ant.raw$status[prb.ant.raw$ENHANCER] <- "enhancer"
cat("  ..promoters\n")
prb.ant.raw$status[grepl("TSS200|5'UTR", prb.ant.raw$UCSC_REFGENE_GROUP)] <- "promoter"

cat("Filtering annotation file by status and selecting relevant columns..\n")
prb.ant <- prb.ant.raw[prb.ant.raw$status != 0,c("TargetID", "CHR", "MAPINFO", "UCSC_REFGENE_NAME", "UCSC_REFGENE_GROUP", "ENHANCER", "status")]

cat("Splitting genes and group...\n")
prb.ant$genes       <- I(strsplit(prb.ant$UCSC_REFGENE_NAME, ";"))
prb.ant$group       <- I(strsplit(prb.ant$UCSC_REFGENE_GROUP, ";"))
cat("Finding and counting unique genes...\n")
prb.ant$unq.genes   <- lapply(prb.ant$genes, unique)
prb.ant$n.unq.genes <- unlist(lapply(prb.ant$unq.genes, length))

# For loop part
cat("Running loop to link probes to genes\n")
probe.loop  <- character(0)
genes.loop  <- character(0)
status.loop <- character(0)

l <- 0
#n <- dim(probe.gene.info)[1]
n <- 50
pb <- txtProgressBar(min=1, max=n, style=3)
for(i in 1:n){
  #setTxtProgressBar(pb, i)
  if(prb.ant$n.unq.genes[i] == 1){
    probe.loop <- c(probe.loop, prb.ant$TargetID[i])
    genes.loop <- c(genes.loop, prb.ant$unq.genes[[i]])
    status.loop <- c(status.loop, prb.ant$status[i])
    l <- l+1
  } else if(prb.ant$n.unq.genes[i] > 1){
    probe  <- prb.ant$TargetID[i]
    genes  <- prb.ant$genes[[i]]
    group  <- prb.ant$group[[i]]
    status <- prb.ant$status[i]

    if(status == "promoter"){
      # Check which genes are promoters
      cat("i =", i, "genes:", genes, "\n")
      cat("i =", i, "group:", group, "\n")
      genes <- unique(genes[grep("TSS200|5'UTR", group)])
      cat("i =", i, "genes:", genes, "\n")
      for(k in 1:length(genes)){
        probe.loop <- c(probe.loop, probe)
        genes.loop <- c(genes.loop, genes[k])
        status.loop <- c(status.loop, status)
        l <- l+1
      }
    } else if (status == "enhancer"){
      # This should be reworked to find nearest genes
      cat("i =", i, "genes:", genes, "\n")
      genes <- unique(genes)
      cat("i =", i, "genes:", genes, "\n")
      for(k in 1:length(genes)){
        probe.loop <- c(probe.loop, probe)
        genes.loop <- c(genes.loop, genes[k])
        status.loop <- c(status.loop, status)
        l <- l+1
      }
    } else {
      # This can be written in a better way, just add unique genes and
      # rep(probe, k) and rep(status, k) to the loop lists
      cat("i =", i, "genes:", genes, "\n")
      genes <- unique(genes)
      cat("i =", i, "genes:", genes, "\n")
      for(k in 1:length(genes)){
        probe.loop <- c(probe.loop, probe)
        genes.loop <- c(genes.loop, genes[k])
        status.loop <- c(status.loop, status)
        l <- l+1
      }
    }
  }
}
