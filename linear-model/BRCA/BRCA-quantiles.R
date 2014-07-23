#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")
load("../Rdata/BRCA/data/BRCA-NMP.Rdata")
load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
load("../Rdata/BRCA/data/BRCA-NEA.Rdata")
load("../Rdata/BRCA/data/BRCA-CEA.Rdata")


# Resize matrices and use only intersecting samples
normal.samples <- intersect(rownames(BRCA.NMP), rownames(BRCA.NEA))
cancer.samples <- intersect(rownames(BRCA.CMP), rownames(BRCA.CEA))

BRCA.NMP.medians        <- colMedians(as.matrix(BRCA.NMP[normal.samples,]), na.rm=T)
names(BRCA.NMP.medians) <- names(BRCA.NMP)
BRCA.CMP.99quant        <- colQuantiles(as.matrix(BRCA.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(BRCA.CMP.99quant) <- names(BRCA.CMP)

BRCA.quantiles <- data.frame(NMP.medians = BRCA.NMP.medians[as.character(BRCA.linked.probes.genes$probes)],
                             CMP.99quant = BRCA.CMP.99quant[as.character(BRCA.linked.probes.genes$probes)])
rownames(BRCA.quantiles)<- paste(BRCA.linked.probes.genes$probe,
                                 BRCA.linked.probes.genes$genes, sep=".")

save(BRCA.quantiles, file="../Rdata/BRCA/calc/BRCA-quantiles.Rdata")
quit(save="no")
