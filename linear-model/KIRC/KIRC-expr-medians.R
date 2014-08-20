#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/KIRC/info/KIRC-linked-probes-genes.Rdata")
load("../Rdata/KIRC/data/KIRC-CMP.Rdata")
load("../Rdata/KIRC/data/KIRC-CEA.Rdata")
if(file.exists("../Rdata/KIRC/data/KIRC-NMP.Rdata") && file.exists("../Rdata/KIRC/data/KIRC-NEA.Rdata")){
  load("../Rdata/KIRC/data/KIRC-NMP.Rdata")
  load("../Rdata/KIRC/data/KIRC-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(KIRC.CMP), rownames(KIRC.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(KIRC.NMP), rownames(KIRC.NEA))
}

chunk.size <- 1000
total.size <- dim(KIRC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

KIRC.expr.nonm   <- rep(NA, total.size)
KIRC.freq.nonm   <- rep(0,  total.size)
KIRC.expr.meth   <- rep(NA, total.size)
KIRC.freq.meth   <- rep(0,  total.size)
KIRC.expr.normal <- rep(NA, total.size)
KIRC.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  KIRC.expr.nonm[index] <- tmp
  KIRC.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  KIRC.expr.meth[index] <- tmp
  KIRC.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- KIRC.NEA[normal.samples, as.character(KIRC.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  KIRC.expr.normal[index] <- tmp
  KIRC.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


KIRC.expr.medians <- data.frame(freq.normal = KIRC.freq.normal, 
                                expr.normal = KIRC.expr.normal,
                                freq.nonm   = KIRC.freq.nonm,
                                expr.nonm   = KIRC.expr.nonm,
                                freq.meth   = KIRC.freq.meth,
                                expr.meth   = KIRC.expr.meth)
rownames(KIRC.expr.medians) <- paste(KIRC.linked.probes.genes$probe, KIRC.linked.probes.genes$genes, sep=".")
save(KIRC.expr.medians, file="../Rdata/KIRC/calc/KIRC-expr-medians.Rdata")
quit(save="no")
