#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LUAD/info/LUAD-linked-probes-genes.Rdata")
load("../Rdata/LUAD/data/LUAD-NMP.Rdata")
load("../Rdata/LUAD/data/LUAD-CMP.Rdata")
load("../Rdata/LUAD/data/LUAD-NEA.Rdata")
load("../Rdata/LUAD/data/LUAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
normal.samples <- intersect(rownames(LUAD.NMP), rownames(LUAD.NEA))
cancer.samples <- intersect(rownames(LUAD.CMP), rownames(LUAD.CEA))

LUAD.NMP.medians        <- colMedians(as.matrix(LUAD.NMP[normal.samples,]), na.rm=T)
names(LUAD.NMP.medians) <- names(LUAD.NMP)
LUAD.CMP.99quant        <- colQuantiles(as.matrix(LUAD.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(LUAD.CMP.99quant) <- names(LUAD.CMP)

LUAD.quantiles <- data.frame(NMP.medians = LUAD.NMP.medians[as.character(LUAD.linked.probes.genes$probes)],
                             CMP.99quant = LUAD.CMP.99quant[as.character(LUAD.linked.probes.genes$probes)])
rownames(LUAD.quantiles)<- paste(LUAD.linked.probes.genes$probe,
                                 LUAD.linked.probes.genes$genes, sep=".")

save(LUAD.quantiles, file="../Rdata/LUAD/calc/LUAD-quantiles.Rdata")
quit(save="no")
