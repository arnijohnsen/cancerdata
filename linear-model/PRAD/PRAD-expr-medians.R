#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/PRAD/info/PRAD-linked-probes-genes.Rdata")
load("../Rdata/PRAD/data/PRAD-CMP.Rdata")
load("../Rdata/PRAD/data/PRAD-CEA.Rdata")
if(file.exists("../Rdata/PRAD/data/PRAD-NMP.Rdata") && file.exists("../Rdata/PRAD/data/PRAD-NEA.Rdata")){
  load("../Rdata/PRAD/data/PRAD-NMP.Rdata")
  load("../Rdata/PRAD/data/PRAD-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(PRAD.CMP), rownames(PRAD.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(PRAD.NMP), rownames(PRAD.NEA))
}

chunk.size <- 1000
total.size <- dim(PRAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

PRAD.expr.nonm   <- rep(NA, total.size)
PRAD.freq.nonm   <- rep(0,  total.size)
PRAD.expr.meth   <- rep(NA, total.size)
PRAD.freq.meth   <- rep(0,  total.size)
PRAD.expr.normal <- rep(NA, total.size)
PRAD.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  PRAD.expr.nonm[index] <- tmp
  PRAD.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  PRAD.expr.meth[index] <- tmp
  PRAD.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- PRAD.NEA[normal.samples, as.character(PRAD.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  PRAD.expr.normal[index] <- tmp
  PRAD.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


PRAD.expr.medians <- data.frame(freq.normal = PRAD.freq.normal, 
                                expr.normal = PRAD.expr.normal,
                                freq.nonm   = PRAD.freq.nonm,
                                expr.nonm   = PRAD.expr.nonm,
                                freq.meth   = PRAD.freq.meth,
                                expr.meth   = PRAD.expr.meth)
rownames(PRAD.expr.medians) <- paste(PRAD.linked.probes.genes$probe, PRAD.linked.probes.genes$genes, sep=".")
save(PRAD.expr.medians, file="../Rdata/PRAD/calc/PRAD-expr-medians.Rdata")
quit(save="no")
