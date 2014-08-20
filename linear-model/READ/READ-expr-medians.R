#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/READ/info/READ-linked-probes-genes.Rdata")
load("../Rdata/READ/data/READ-CMP.Rdata")
load("../Rdata/READ/data/READ-CEA.Rdata")
if(file.exists("../Rdata/READ/data/READ-NMP.Rdata") && file.exists("../Rdata/READ/data/READ-NEA.Rdata")){
  load("../Rdata/READ/data/READ-NMP.Rdata")
  load("../Rdata/READ/data/READ-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(READ.CMP), rownames(READ.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(READ.NMP), rownames(READ.NEA))
}

chunk.size <- 1000
total.size <- dim(READ.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

READ.expr.nonm   <- rep(NA, total.size)
READ.freq.nonm   <- rep(0,  total.size)
READ.expr.meth   <- rep(NA, total.size)
READ.freq.meth   <- rep(0,  total.size)
READ.expr.normal <- rep(NA, total.size)
READ.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  READ.expr.nonm[index] <- tmp
  READ.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  READ.expr.meth[index] <- tmp
  READ.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- READ.NEA[normal.samples, as.character(READ.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  READ.expr.normal[index] <- tmp
  READ.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


READ.expr.medians <- data.frame(freq.normal = READ.freq.normal, 
                                expr.normal = READ.expr.normal,
                                freq.nonm   = READ.freq.nonm,
                                expr.nonm   = READ.expr.nonm,
                                freq.meth   = READ.freq.meth,
                                expr.meth   = READ.expr.meth)
rownames(READ.expr.medians) <- paste(READ.linked.probes.genes$probe, READ.linked.probes.genes$genes, sep=".")
save(READ.expr.medians, file="../Rdata/READ/calc/READ-expr-medians.Rdata")
quit(save="no")
