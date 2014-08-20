#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/OV/info/OV-linked-probes-genes.Rdata")
load("../Rdata/OV/data/OV-CMP.Rdata")
load("../Rdata/OV/data/OV-CEA.Rdata")
if(file.exists("../Rdata/OV/data/OV-NMP.Rdata") && file.exists("../Rdata/OV/data/OV-NEA.Rdata")){
  load("../Rdata/OV/data/OV-NMP.Rdata")
  load("../Rdata/OV/data/OV-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(OV.CMP), rownames(OV.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(OV.NMP), rownames(OV.NEA))
}

chunk.size <- 1000
total.size <- dim(OV.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

OV.expr.nonm   <- rep(NA, total.size)
OV.freq.nonm   <- rep(0,  total.size)
OV.expr.meth   <- rep(NA, total.size)
OV.freq.meth   <- rep(0,  total.size)
OV.expr.normal <- rep(NA, total.size)
OV.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  OV.expr.nonm[index] <- tmp
  OV.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  OV.expr.meth[index] <- tmp
  OV.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- OV.NEA[normal.samples, as.character(OV.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  OV.expr.normal[index] <- tmp
  OV.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


OV.expr.medians <- data.frame(freq.normal = OV.freq.normal, 
                                expr.normal = OV.expr.normal,
                                freq.nonm   = OV.freq.nonm,
                                expr.nonm   = OV.expr.nonm,
                                freq.meth   = OV.freq.meth,
                                expr.meth   = OV.expr.meth)
rownames(OV.expr.medians) <- paste(OV.linked.probes.genes$probe, OV.linked.probes.genes$genes, sep=".")
save(OV.expr.medians, file="../Rdata/OV/calc/OV-expr-medians.Rdata")
quit(save="no")
