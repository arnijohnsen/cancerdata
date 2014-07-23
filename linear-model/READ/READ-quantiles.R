#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/READ/info/READ-linked-probes-genes.Rdata")
load("../Rdata/READ/data/READ-NMP.Rdata")
load("../Rdata/READ/data/READ-CMP.Rdata")
load("../Rdata/READ/data/READ-NEA.Rdata")
load("../Rdata/READ/data/READ-CEA.Rdata")


# Resize matrices and use only intersecting samples
normal.samples <- intersect(rownames(READ.NMP), rownames(READ.NEA))
cancer.samples <- intersect(rownames(READ.CMP), rownames(READ.CEA))

READ.NMP.medians        <- colMedians(as.matrix(READ.NMP[normal.samples,]), na.rm=T)
names(READ.NMP.medians) <- names(READ.NMP)
READ.CMP.99quant        <- colQuantiles(as.matrix(READ.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(READ.CMP.99quant) <- names(READ.CMP)

READ.quantiles <- data.frame(NMP.medians = READ.NMP.medians[as.character(READ.linked.probes.genes$probes)],
                             CMP.99quant = READ.CMP.99quant[as.character(READ.linked.probes.genes$probes)])
rownames(READ.quantiles)<- paste(READ.linked.probes.genes$probe,
                                 READ.linked.probes.genes$genes, sep=".")

save(READ.quantiles, file="../Rdata/READ/calc/READ-quantiles.Rdata")
quit(save="no")
