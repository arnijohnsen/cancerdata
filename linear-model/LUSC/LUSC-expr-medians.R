#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LUSC/info/LUSC-linked-probes-genes.Rdata")
load("../Rdata/LUSC/data/LUSC-CMP.Rdata")
load("../Rdata/LUSC/data/LUSC-CEA.Rdata")
if(file.exists("../Rdata/LUSC/data/LUSC-NMP.Rdata") && file.exists("../Rdata/LUSC/data/LUSC-NEA.Rdata")){
  load("../Rdata/LUSC/data/LUSC-NMP.Rdata")
  load("../Rdata/LUSC/data/LUSC-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LUSC.CMP), rownames(LUSC.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(LUSC.NMP), rownames(LUSC.NEA))
}

chunk.size <- 1000
total.size <- dim(LUSC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LUSC.expr.nonm   <- rep(NA, total.size)
LUSC.freq.nonm   <- rep(0,  total.size)
LUSC.expr.meth   <- rep(NA, total.size)
LUSC.freq.meth   <- rep(0,  total.size)
LUSC.expr.normal <- rep(NA, total.size)
LUSC.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- LUSC.CMP[cancer.samples, as.character(LUSC.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- LUSC.CEA[cancer.samples, as.character(LUSC.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LUSC.expr.nonm[index] <- tmp
  LUSC.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- LUSC.CMP[cancer.samples, as.character(LUSC.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- LUSC.CEA[cancer.samples, as.character(LUSC.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LUSC.expr.meth[index] <- tmp
  LUSC.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- LUSC.NEA[normal.samples, as.character(LUSC.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LUSC.expr.normal[index] <- tmp
  LUSC.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


LUSC.expr.medians <- data.frame(freq.normal = LUSC.freq.normal, 
                                expr.normal = LUSC.expr.normal,
                                freq.nonm   = LUSC.freq.nonm,
                                expr.nonm   = LUSC.expr.nonm,
                                freq.meth   = LUSC.freq.meth,
                                expr.meth   = LUSC.expr.meth)
rownames(LUSC.expr.medians) <- paste(LUSC.linked.probes.genes$probe, LUSC.linked.probes.genes$genes, sep=".")
save(LUSC.expr.medians, file="../Rdata/LUSC/calc/LUSC-expr-medians.Rdata")
quit(save="no")
