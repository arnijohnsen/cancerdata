#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LIHC/info/LIHC-linked-probes-genes.Rdata")
load("../Rdata/LIHC/data/LIHC-CMP.Rdata")
load("../Rdata/LIHC/data/LIHC-CEA.Rdata")
if(file.exists("../Rdata/LIHC/data/LIHC-NMP.Rdata") && file.exists("../Rdata/LIHC/data/LIHC-NEA.Rdata")){
  load("../Rdata/LIHC/data/LIHC-NMP.Rdata")
  load("../Rdata/LIHC/data/LIHC-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LIHC.CMP), rownames(LIHC.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(LIHC.NMP), rownames(LIHC.NEA))
}

chunk.size <- 1000
total.size <- dim(LIHC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LIHC.expr.nonm   <- rep(NA, total.size)
LIHC.freq.nonm   <- rep(0,  total.size)
LIHC.expr.meth   <- rep(NA, total.size)
LIHC.freq.meth   <- rep(0,  total.size)
LIHC.expr.normal <- rep(NA, total.size)
LIHC.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LIHC.expr.nonm[index] <- tmp
  LIHC.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LIHC.expr.meth[index] <- tmp
  LIHC.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- LIHC.NEA[normal.samples, as.character(LIHC.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LIHC.expr.normal[index] <- tmp
  LIHC.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


LIHC.expr.medians <- data.frame(freq.normal = LIHC.freq.normal, 
                                expr.normal = LIHC.expr.normal,
                                freq.nonm   = LIHC.freq.nonm,
                                expr.nonm   = LIHC.expr.nonm,
                                freq.meth   = LIHC.freq.meth,
                                expr.meth   = LIHC.expr.meth)
rownames(LIHC.expr.medians) <- paste(LIHC.linked.probes.genes$probe, LIHC.linked.probes.genes$genes, sep=".")
save(LIHC.expr.medians, file="../Rdata/LIHC/calc/LIHC-expr-medians.Rdata")
quit(save="no")
