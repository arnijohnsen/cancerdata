#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/COAD/info/COAD-linked-probes-genes.Rdata")
load("../Rdata/COAD/data/COAD-CMP.Rdata")
load("../Rdata/COAD/data/COAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(COAD.CMP), rownames(COAD.CEA))

COAD.CMP.99quant        <- colQuantiles(as.matrix(COAD.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(COAD.CMP.99quant) <- names(COAD.CMP)

COAD.quantiles <- data.frame(CMP.99quant = COAD.CMP.99quant[as.character(COAD.linked.probes.genes$probes)])
rownames(COAD.quantiles)<- paste(COAD.linked.probes.genes$probe,
                                 COAD.linked.probes.genes$genes, sep=".")

save(COAD.quantiles, file="../Rdata/COAD/calc/COAD-quantiles.Rdata")
quit(save="no")
