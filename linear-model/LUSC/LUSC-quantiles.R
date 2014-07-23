#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LUSC/info/LUSC-linked-probes-genes.Rdata")
load("../Rdata/LUSC/data/LUSC-NMP.Rdata")
load("../Rdata/LUSC/data/LUSC-CMP.Rdata")
load("../Rdata/LUSC/data/LUSC-NEA.Rdata")
load("../Rdata/LUSC/data/LUSC-CEA.Rdata")


# Resize matrices and use only intersecting samples
normal.samples <- intersect(rownames(LUSC.NMP), rownames(LUSC.NEA))
cancer.samples <- intersect(rownames(LUSC.CMP), rownames(LUSC.CEA))

LUSC.NMP.medians        <- colMedians(as.matrix(LUSC.NMP[normal.samples,]), na.rm=T)
names(LUSC.NMP.medians) <- names(LUSC.NMP)
LUSC.CMP.99quant        <- colQuantiles(as.matrix(LUSC.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(LUSC.CMP.99quant) <- names(LUSC.CMP)

LUSC.quantiles <- data.frame(NMP.medians = LUSC.NMP.medians[as.character(LUSC.linked.probes.genes$probes)],
                             CMP.99quant = LUSC.CMP.99quant[as.character(LUSC.linked.probes.genes$probes)])
rownames(LUSC.quantiles)<- paste(LUSC.linked.probes.genes$probe,
                                 LUSC.linked.probes.genes$genes, sep=".")

save(LUSC.quantiles, file="../Rdata/LUSC/calc/LUSC-quantiles.Rdata")
quit(save="no")
