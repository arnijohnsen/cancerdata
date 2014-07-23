#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/READ/info/READ-linked-probes-genes.Rdata")
load("../Rdata/READ/data/READ-CMP.Rdata")
load("../Rdata/READ/data/READ-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(READ.CMP), rownames(READ.CEA))

chunk.size <- 1000
total.size <- dim(READ.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

READ.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  READ.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(READ.nm.expr) <- paste(READ.linked.probes.genes$probe,
                             READ.linked.probes.genes$genes, sep=".")
save(READ.nm.expr, file="../Rdata/READ/calc/READ-nm-expr.Rdata")
quit(save="no")
