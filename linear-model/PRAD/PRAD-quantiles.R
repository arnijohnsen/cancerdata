#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/PRAD/info/PRAD-linked-probes-genes.Rdata")
load("../Rdata/PRAD/data/PRAD-NMP.Rdata")
load("../Rdata/PRAD/data/PRAD-CMP.Rdata")
load("../Rdata/PRAD/data/PRAD-NEA.Rdata")
load("../Rdata/PRAD/data/PRAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
normal.samples <- intersect(rownames(PRAD.NMP), rownames(PRAD.NEA))
cancer.samples <- intersect(rownames(PRAD.CMP), rownames(PRAD.CEA))

PRAD.NMP.medians        <- colMedians(as.matrix(PRAD.NMP[normal.samples,]), na.rm=T)
names(PRAD.NMP.medians) <- names(PRAD.NMP)
PRAD.CMP.99quant        <- colQuantiles(as.matrix(PRAD.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(PRAD.CMP.99quant) <- names(PRAD.CMP)

PRAD.quantiles <- data.frame(NMP.medians = PRAD.NMP.medians[as.character(PRAD.linked.probes.genes$probes)],
                             CMP.99quant = PRAD.CMP.99quant[as.character(PRAD.linked.probes.genes$probes)])
rownames(PRAD.quantiles)<- paste(PRAD.linked.probes.genes$probe,
                                 PRAD.linked.probes.genes$genes, sep=".")

save(PRAD.quantiles, file="../Rdata/PRAD/calc/PRAD-quantiles.Rdata")
quit(save="no")
