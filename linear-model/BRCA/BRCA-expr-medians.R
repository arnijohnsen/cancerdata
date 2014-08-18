#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")
load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
load("../Rdata/BRCA/data/BRCA-CEA.Rdata")
load("../Rdata/BRCA/data/BRCA-NMP.Rdata")
load("../Rdata/BRCA/data/BRCA-NEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(BRCA.CMP), rownames(BRCA.CEA))
normal.samples <- intersect(rownames(BRCA.NMP), rownames(BRCA.NEA))

chunk.size <- 1000
total.size <- dim(BRCA.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

BRCA.expr.nonm   <- rep(0, total.size)
BRCA.freq.nonm   <- rep(0, total.size)
BRCA.expr.meth   <- rep(0, total.size)
BRCA.freq.meth   <- rep(0, total.size)
BRCA.expr.normal <- rep(0, total.size)
BRCA.freq.normal <- rep(0, total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  BRCA.expr.nonm[index] <- tmp
  BRCA.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  BRCA.expr.meth[index] <- tmp
  BRCA.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- BRCA.NEA[normal.samples, as.character(BRCA.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  BRCA.expr.normal[index] <- tmp
  BRCA.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

BRCA.expr.medians <- data.frame(freq.normal = BRCA.freq.normal, 
                                expr.normal = BRCA.expr.normal,
                                freq.nonm   = BRCA.freq.nonm,
                                expr.nonm   = BRCA.expr.nonm,
                                freq.meth   = BRCA.freq.meth,
                                expr.meth   = BRCA.expr.meth)
rownames(BRCA.expr.medians) <- paste(BRCA.linked.probes.genes$probe, BRCA.linked.probes.genes$genes, sep=".")
save(BRCA.expr.medians, file="../Rdata/BRCA/calc/BRCA-expr-medians.Rdata")
quit(save="no")
