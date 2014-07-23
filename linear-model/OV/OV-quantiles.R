#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/OV/info/OV-linked-probes-genes.Rdata")
load("../Rdata/OV/data/OV-CMP.Rdata")
load("../Rdata/OV/data/OV-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(OV.CMP), rownames(OV.CEA))

OV.CMP.99quant        <- colQuantiles(as.matrix(OV.CMP[cancer.samples,]), probs=0.99, na.rm=T)
names(OV.CMP.99quant) <- names(OV.CMP)

OV.quantiles <- data.frame(CMP.99quant = OV.CMP.99quant[as.character(OV.linked.probes.genes$probes)])
rownames(OV.quantiles)<- paste(OV.linked.probes.genes$probe,
                                 OV.linked.probes.genes$genes, sep=".")

save(OV.quantiles, file="../Rdata/OV/calc/OV-quantiles.Rdata")
quit(save="no")
