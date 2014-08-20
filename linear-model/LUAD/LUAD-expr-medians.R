#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LUAD/info/LUAD-linked-probes-genes.Rdata")
load("../Rdata/LUAD/data/LUAD-CMP.Rdata")
load("../Rdata/LUAD/data/LUAD-CEA.Rdata")
if(file.exists("../Rdata/LUAD/data/LUAD-NMP.Rdata") && file.exists("../Rdata/LUAD/data/LUAD-NEA.Rdata")){
  load("../Rdata/LUAD/data/LUAD-NMP.Rdata")
  load("../Rdata/LUAD/data/LUAD-NEA.Rdata")
  has.normal <- TRUE
} else {
  has.normal <- FALSE
}


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LUAD.CMP), rownames(LUAD.CEA))
if(has.normal){
  normal.samples <- intersect(rownames(LUAD.NMP), rownames(LUAD.NEA))
}

chunk.size <- 1000
total.size <- dim(LUAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LUAD.expr.nonm   <- rep(NA, total.size)
LUAD.freq.nonm   <- rep(0,  total.size)
LUAD.expr.meth   <- rep(NA, total.size)
LUAD.freq.meth   <- rep(0,  total.size)
LUAD.expr.normal <- rep(NA, total.size)
LUAD.freq.normal <- rep(0,  total.size)

cat("Running loop for non methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] <= (1/3)
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LUAD.expr.nonm[index] <- tmp
  LUAD.freq.nonm[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

cat("Running loop for methylated\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] > (1/3)
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LUAD.expr.meth[index] <- tmp
  LUAD.freq.meth[index] <- apply(x, 2, function(x) {sum(!is.na(x))})
}
cat("\n")

if(has.normal){
cat("Running loop for normal\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))

  y <- LUAD.NEA[normal.samples, as.character(LUAD.linked.probes.genes$genes[index])]
  tmp <- colMedians(as.matrix(y), na.rm=T)
  tmp[is.nan(tmp)] <- NA
  LUAD.expr.normal[index] <- tmp
  LUAD.freq.normal[index] <- apply(y, 2, function(x) {sum(!is.na(x))})
}
cat("\n")
}


LUAD.expr.medians <- data.frame(freq.normal = LUAD.freq.normal, 
                                expr.normal = LUAD.expr.normal,
                                freq.nonm   = LUAD.freq.nonm,
                                expr.nonm   = LUAD.expr.nonm,
                                freq.meth   = LUAD.freq.meth,
                                expr.meth   = LUAD.expr.meth)
rownames(LUAD.expr.medians) <- paste(LUAD.linked.probes.genes$probe, LUAD.linked.probes.genes$genes, sep=".")
save(LUAD.expr.medians, file="../Rdata/LUAD/calc/LUAD-expr-medians.Rdata")
quit(save="no")
