#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LIHC/info/LIHC-linked-probes-genes.Rdata")
load("../Rdata/LIHC/data/LIHC-NMP.Rdata")
load("../Rdata/LIHC/data/LIHC-CMP.Rdata")
load("../Rdata/LIHC/data/LIHC-NEA.Rdata")
load("../Rdata/LIHC/data/LIHC-CEA.Rdata")

# Resize matrices and use only intersecting samples
normal.samples <- intersect(rownames(LIHC.NMP), rownames(LIHC.NEA))
cancer.samples <- intersect(rownames(LIHC.CMP), rownames(LIHC.CEA))

LIHC.NMP.medians        <- colMedians(as.matrix(LIHC.NMP[normal.samples,]), na.rm=T)
names(LIHC.NMP.medians) <- names(LIHC.NMP)
LIHC.CMP.99quant        <- colQuantiles(as.matrix(LIHC.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(LIHC.CMP.99quant) <- names(LIHC.CMP)

LIHC.quantiles <- data.frame(NMP.medians = LIHC.NMP.medians[as.character(LIHC.linked.probes.genes$probes)],
                             CMP.99quant = LIHC.CMP.99quant[as.character(LIHC.linked.probes.genes$probes)])
rownames(LIHC.quantiles)<- paste(LIHC.linked.probes.genes$probe,
                                 LIHC.linked.probes.genes$genes, sep=".")

save(LIHC.quantiles, file="../Rdata/LIHC/calc/LIHC-quantiles.Rdata")
quit(save="no")
