#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/KIRC/info/KIRC-linked-probes-genes.Rdata")
load("../Rdata/KIRC/data/KIRC-CMP.Rdata")
load("../Rdata/KIRC/data/KIRC-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(KIRC.CMP), rownames(KIRC.CEA))

KIRC.CMP.99quant        <- colQuantiles(as.matrix(KIRC.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(KIRC.CMP.99quant) <- names(KIRC.CMP)

KIRC.quantiles <- data.frame(CMP.99quant = KIRC.CMP.99quant[as.character(KIRC.linked.probes.genes$probes)])
rownames(KIRC.quantiles)<- paste(KIRC.linked.probes.genes$probe,
                                 KIRC.linked.probes.genes$genes, sep=".")

save(KIRC.quantiles, file="../Rdata/KIRC/calc/KIRC-quantiles.Rdata")
quit(save="no")
