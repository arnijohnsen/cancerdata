#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LUSC/info/LUSC-linked-probes-genes.Rdata")
load("../Rdata/LUSC/data/LUSC-CMP.Rdata")
load("../Rdata/LUSC/data/LUSC-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LUSC.CMP), rownames(LUSC.CEA))

chunk.size <- 1000
total.size <- dim(LUSC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LUSC.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LUSC.CMP[cancer.samples, as.character(LUSC.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- LUSC.CEA[cancer.samples, as.character(LUSC.linked.probes.genes$genes[index])]*x
  LUSC.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(LUSC.nm.expr) <- paste(LUSC.linked.probes.genes$probe,
                             LUSC.linked.probes.genes$genes, sep=".")
save(LUSC.nm.expr, file="../Rdata/LUSC/calc/LUSC-nm-expr.Rdata")
quit(save="no")
