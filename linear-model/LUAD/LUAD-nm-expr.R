#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LUAD/info/LUAD-linked-probes-genes.Rdata")
load("../Rdata/LUAD/data/LUAD-CMP.Rdata")
load("../Rdata/LUAD/data/LUAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LUAD.CMP), rownames(LUAD.CEA))

chunk.size <- 1000
total.size <- dim(LUAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LUAD.nm.expr <- rep(0, total.size)

cat("Running loop\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] < 0.2
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  LUAD.nm.expr[index] <-  colQuantiles(y, probs=0.10, na.rm=T)
}
cat("\n")

names(LUAD.nm.expr) <- paste(LUAD.linked.probes.genes$probe,
                             LUAD.linked.probes.genes$genes, sep=".")
save(LUAD.nm.expr, file="../Rdata/LUAD/calc/LUAD-nm-expr.Rdata")
quit(save="no")
