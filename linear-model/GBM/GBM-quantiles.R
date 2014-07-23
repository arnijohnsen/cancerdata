#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/GBM/info/GBM-linked-probes-genes.Rdata")
load("../Rdata/GBM/data/GBM-CMP.Rdata")
load("../Rdata/GBM/data/GBM-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(GBM.CMP), rownames(GBM.CEA))

GBM.CMP.99quant        <- colQuantiles(as.matrix(GBM.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(GBM.CMP.99quant) <- names(GBM.CMP)

GBM.quantiles <- data.frame(CMP.99quant = GBM.CMP.99quant[as.character(GBM.linked.probes.genes$probes)])
rownames(GBM.quantiles)<- paste(GBM.linked.probes.genes$probe,
                                 GBM.linked.probes.genes$genes, sep=".")

save(GBM.quantiles, file="../Rdata/GBM/calc/GBM-quantiles.Rdata")
quit(save="no")
