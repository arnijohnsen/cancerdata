cat("Reading probes annotation file..\n")
probe.annotations <- read.table("../rawdata/AnnotationFiles/GenomeStudioProbeAnnotations.txt", header=TRUE, sep="\t", quote="\"", stringsAsFactors=FALSE)

cat("Setting probe status..\n")
probe.annotations$status <- 0

cat("  ..gene bodies\n")
is.body   <- !grepl("TSS200|TSS1500|1stExon|5'UTR|3'UTR|^$", probe.annotations$UCSC_REFGENE_GROUP)
is.island <-  grepl("Island", probe.annotations$RELATION_TO_UCSC_CPG_ISLAND)
is.shore  <-  grepl("Shore",  probe.annotations$RELATION_TO_UCSC_CPG_ISLAND)
is.none   <- !is.island & !is.shore
probe.annotations$status[is.body & is.island] <- "body.island"
probe.annotations$status[is.body & is.shore]  <- "body.shore"
probe.annotations$status[is.body & is.none]   <- "body.none"

cat("  ..enhancers\n")
probe.annotations$status[probe.annotations$ENHANCER] <- "enhancer"
cat("  ..promoters\n")
probe.annotations$status[grepl("TSS200|5'UTR", probe.annotations$UCSC_REFGENE_GROUP)] <- "promoter"

chromosomes <- c(1:22, "X", "Y")
body.island.probes <- as.list(chromosomes)
body.shore.probes  <- as.list(chromosomes)
body.none.probes   <- as.list(chromosomes)
promoter.probes    <- as.list(chromosomes)
enhancer.probes    <- as.list(chromosomes)
names(body.island.probes) <- paste("chr.",chromosomes, sep="")
names(body.shore.probes ) <- paste("chr.",chromosomes, sep="")
names(body.none.probes  ) <- paste("chr.",chromosomes, sep="")
names(promoter.probes   ) <- paste("chr.",chromosomes, sep="")
names(enhancer.probes   ) <- paste("chr.",chromosomes, sep="")
pb <- txtProgressBar(min=1, max=length(chromosomes), style=3)
for(i in 1:length(chromosomes)){
  setTxtProgressBar(pb,i)
  body.island.probes[[paste("chr.",chromosomes[i],sep="")]] <- probe.annotations$TargetID[probe.annotations$CHR == chromosomes[i] & probe.annotations$status == "body.island"]
  body.shore.probes[[ paste("chr.",chromosomes[i],sep="")]] <- probe.annotations$TargetID[probe.annotations$CHR == chromosomes[i] & probe.annotations$status == "body.shore" ]
  body.none.probes[[  paste("chr.",chromosomes[i],sep="")]] <- probe.annotations$TargetID[probe.annotations$CHR == chromosomes[i] & probe.annotations$status == "body.none"  ]
  promoter.probes[[   paste("chr.",chromosomes[i],sep="")]] <- probe.annotations$TargetID[probe.annotations$CHR == chromosomes[i] & probe.annotations$status == "promoter"   ]
  enhancer.probes[[   paste("chr.",chromosomes[i],sep="")]] <- probe.annotations$TargetID[probe.annotations$CHR == chromosomes[i] & probe.annotations$status == "enhancer"   ]
}
cat("\n")
save(body.island.probes, file="../Rdata/BRCA-new/methylation/probe-annotations/body-island-probes.Rdata")
save(body.shore.probes , file="../Rdata/BRCA-new/methylation/probe-annotations/body-shore-probes.Rdata")
save(body.none.probes  , file="../Rdata/BRCA-new/methylation/probe-annotations/body-none-probes.Rdata")
save(promoter.probes   , file="../Rdata/BRCA-new/methylation/probe-annotations/promoter-probes.Rdata")
save(enhancer.probes   , file="../Rdata/BRCA-new/methylation/probe-annotations/enhancer-probes.Rdata")
quit(save="no")
