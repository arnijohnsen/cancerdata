#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LIHC/info/LIHC-linked-probes-genes.Rdata")
load("../Rdata/LIHC/data/LIHC-CMP.Rdata")
load("../Rdata/LIHC/data/LIHC-CEA.Rdata")

LIHC.linked.probes.genes <- LIHC.linked.probes.genes[-53182,]

# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LIHC.CMP), rownames(LIHC.CEA))

chunk.size <- 1000
total.size <- dim(LIHC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LIHC.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  LIHC.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(LIHC.nm.expr) <- paste(LIHC.linked.probes.genes$probe,
                             LIHC.linked.probes.genes$genes, sep=".")
save(LIHC.nm.expr, file="../Rdata/LIHC/calc/LIHC-nm-expr.Rdata")
quit(save="no")
