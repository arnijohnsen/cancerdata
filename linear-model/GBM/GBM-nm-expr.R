#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/GBM/info/GBM-linked-probes-genes.Rdata")
load("../Rdata/GBM/data/GBM-CMP.Rdata")
load("../Rdata/GBM/data/GBM-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(GBM.CMP), rownames(GBM.CEA))

chunk.size <- 1000
total.size <- dim(GBM.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

GBM.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  GBM.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(GBM.nm.expr) <- paste(GBM.linked.probes.genes$probe,
                             GBM.linked.probes.genes$genes, sep=".")
save(GBM.nm.expr, file="../Rdata/GBM/calc/GBM-nm-expr.Rdata")
quit(save="no")
