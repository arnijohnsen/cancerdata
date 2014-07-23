#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/COAD/info/COAD-linked-probes-genes.Rdata")
load("../Rdata/COAD/data/COAD-CMP.Rdata")
load("../Rdata/COAD/data/COAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(COAD.CMP), rownames(COAD.CEA))

chunk.size <- 1000
total.size <- dim(COAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

COAD.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  COAD.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(COAD.nm.expr) <- paste(COAD.linked.probes.genes$probe,
                             COAD.linked.probes.genes$genes, sep=".")
save(COAD.nm.expr, file="../Rdata/COAD/calc/COAD-nm-expr.Rdata")
quit(save="no")
