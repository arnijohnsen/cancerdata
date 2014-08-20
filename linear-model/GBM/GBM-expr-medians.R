#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/GBM/info/GBM-linked-probes-genes.Rdata")
load("../Rdata/GBM/data/GBM-CMP.Rdata")
load("../Rdata/GBM/data/GBM-CEA.Rdata")
if(file.exists("../Rdata/GBM/data/GBM-NMP.Rdata") && file.exists("../Rdata/GBM/data/GBM-NEA.Rdata")){
  load("../Rdata/GBM/data/GBM-NMP.Rdata")
  load("../Rdata/GBM/data/GBM-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(GBM.CMP), rownames(GBM.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(GBM.NMP), rownames(GBM.NEA))
}

chunk.size <- 1000
total.size <- dim(GBM.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

GBM.expr.nonm   <- rep(NA, total.size)
GBM.freq.nonm   <- rep(0,  total.size)
GBM.expr.meth   <- rep(NA, total.size)
GBM.freq.meth   <- rep(0,  total.size)
GBM.expr.normal <- rep(NA, total.size)
GBM.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  GBM.expr.nonm[index] <- tmp
  GBM.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  GBM.expr.meth[index] <- tmp
  GBM.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- GBM.NEA[normal.samples, as.character(GBM.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  GBM.expr.normal[index] <- tmp
  GBM.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


GBM.expr.medians <- data.frame(freq.normal = GBM.freq.normal, 
                                expr.normal = GBM.expr.normal,
                                freq.nonm   = GBM.freq.nonm,
                                expr.nonm   = GBM.expr.nonm,
                                freq.meth   = GBM.freq.meth,
                                expr.meth   = GBM.expr.meth)
rownames(GBM.expr.medians) <- paste(GBM.linked.probes.genes$probe, GBM.linked.probes.genes$genes, sep=".")
save(GBM.expr.medians, file="../Rdata/GBM/calc/GBM-expr-medians.Rdata")
quit(save="no")
