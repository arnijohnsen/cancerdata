#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/PRAD/info/PRAD-linked-probes-genes.Rdata")
load("../Rdata/PRAD/data/PRAD-CMP.Rdata")
load("../Rdata/PRAD/data/PRAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(PRAD.CMP), rownames(PRAD.CEA))

chunk.size <- 1000
total.size <- dim(PRAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

PRAD.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  PRAD.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(PRAD.nm.expr) <- paste(PRAD.linked.probes.genes$probe,
                             PRAD.linked.probes.genes$genes, sep=".")
save(PRAD.nm.expr, file="../Rdata/PRAD/calc/PRAD-nm-expr.Rdata")
quit(save="no")
