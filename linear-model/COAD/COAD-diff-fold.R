#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/COAD/info/COAD-linked-probes-genes.Rdata")
load("../Rdata/COAD/data/COAD-CMP.Rdata")
load("../Rdata/COAD/data/COAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(COAD.CMP), rownames(COAD.CEA))

chunk.size <- 1000
total.size <- dim(COAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

COAD.diff.25 <- rep(0, total.size)
COAD.fold.25 <- rep(0, total.size)
COAD.diff.33 <- rep(0, total.size)
COAD.fold.33 <- rep(0, total.size)
COAD.diff.40 <- rep(0, total.size)
COAD.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  COAD.diff.25[index] <- m2-m1
  COAD.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  COAD.diff.33[index] <- m2-m1
  COAD.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- COAD.CMP[cancer.samples, as.character(COAD.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- COAD.CEA[cancer.samples, as.character(COAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  COAD.diff.40[index] <- m2-m1
  COAD.fold.40[index] <- m2/m1
}
cat("\n")

COAD.diff.fold <- data.frame(diff.25 = COAD.diff.25,
                           diff.33 = COAD.diff.33,
                           diff.40 = COAD.diff.40,
                           fold.25 = COAD.fold.25,
                           fold.33 = COAD.fold.33,
                           fold.40 = COAD.fold.40)
rownames(COAD.diff.fold) <- paste(COAD.linked.probes.genes$probe, COAD.linked.probes.genes$genes, sep=".")
COAD.diff.fold[is.nan(as.matrix(COAD.diff.fold))] <- NA
save(COAD.diff.fold, file="../Rdata/COAD/calc/COAD-diff-fold.Rdata")
quit(save="no")
