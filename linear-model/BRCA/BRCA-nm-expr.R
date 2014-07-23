#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")
load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
load("../Rdata/BRCA/data/BRCA-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(BRCA.CMP), rownames(BRCA.CEA))

chunk.size <- 1000
total.size <- dim(BRCA.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

BRCA.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  BRCA.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(BRCA.nm.expr) <- paste(BRCA.linked.probes.genes$probe,
                             BRCA.linked.probes.genes$genes, sep=".")
save(BRCA.nm.expr, file="../Rdata/BRCA/calc/BRCA-nm-expr.Rdata")
quit(save="no")
