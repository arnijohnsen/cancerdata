#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/KIRC/info/KIRC-linked-probes-genes.Rdata")
load("../Rdata/KIRC/data/KIRC-CMP.Rdata")
load("../Rdata/KIRC/data/KIRC-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(KIRC.CMP), rownames(KIRC.CEA))

chunk.size <- 1000
total.size <- dim(KIRC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

KIRC.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  KIRC.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(KIRC.nm.expr) <- paste(KIRC.linked.probes.genes$probe,
                             KIRC.linked.probes.genes$genes, sep=".")
save(KIRC.nm.expr, file="../Rdata/KIRC/calc/KIRC-nm-expr.Rdata")
quit(save="no")
