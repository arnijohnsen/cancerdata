#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/OV/info/OV-linked-probes-genes.Rdata")
load("../Rdata/OV/data/OV-CMP.Rdata")
load("../Rdata/OV/data/OV-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(OV.CMP), rownames(OV.CEA))

chunk.size <- 1000
total.size <- dim(OV.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

OV.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  OV.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(OV.nm.expr) <- paste(OV.linked.probes.genes$probe,
                             OV.linked.probes.genes$genes, sep=".")
save(OV.nm.expr, file="../Rdata/OV/calc/OV-nm-expr.Rdata")
quit(save="no")
