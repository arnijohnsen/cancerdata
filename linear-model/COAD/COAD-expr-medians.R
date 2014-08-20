#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/COAD/info/COAD-linked-probes-genes.Rdata")
load("../Rdata/COAD/data/COAD-CMP.Rdata")
load("../Rdata/COAD/data/COAD-CEA.Rdata")
if(file.exists("../Rdata/COAD/data/COAD-NMP.Rdata") && file.exists("../Rdata/COAD/data/COAD-NEA.Rdata")){
  load("../Rdata/COAD/data/COAD-NMP.Rdata")
  load("../Rdata/COAD/data/COAD-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(COAD.CMP), rownames(COAD.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(COAD.NMP), rownames(COAD.NEA))
}

chunk.size <- 1000
total.size <- dim(COAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

COAD.expr.nonm   <- rep(NA, total.size)
COAD.freq.nonm   <- rep(0,  total.size)
COAD.expr.meth   <- rep(NA, total.size)
COAD.freq.meth   <- rep(0,  total.size)
COAD.expr.normal <- rep(NA, total.size)
COAD.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  COAD.expr.nonm[index] <- tmp
  COAD.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  COAD.expr.meth[index] <- tmp
  COAD.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- COAD.NEA[normal.samples, as.character(COAD.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  COAD.expr.normal[index] <- tmp
  COAD.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


COAD.expr.medians <- data.frame(freq.normal = COAD.freq.normal, 
                                expr.normal = COAD.expr.normal,
                                freq.nonm   = COAD.freq.nonm,
                                expr.nonm   = COAD.expr.nonm,
                                freq.meth   = COAD.freq.meth,
                                expr.meth   = COAD.expr.meth)
rownames(COAD.expr.medians) <- paste(COAD.linked.probes.genes$probe, COAD.linked.probes.genes$genes, sep=".")
save(COAD.expr.medians, file="../Rdata/COAD/calc/COAD-expr-medians.Rdata")
quit(save="no")
