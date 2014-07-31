#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LIHC/info/LIHC-linked-probes-genes.Rdata")
load("../Rdata/LIHC/data/LIHC-CMP.Rdata")
load("../Rdata/LIHC/data/LIHC-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LIHC.CMP), rownames(LIHC.CEA))

chunk.size <- 1000
total.size <- dim(LIHC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LIHC.diff.25 <- rep(0, total.size)
LIHC.fold.25 <- rep(0, total.size)
LIHC.diff.33 <- rep(0, total.size)
LIHC.fold.33 <- rep(0, total.size)
LIHC.diff.40 <- rep(0, total.size)
LIHC.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  LIHC.diff.25[index] <- m2-m1
  LIHC.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  LIHC.diff.33[index] <- m2-m1
  LIHC.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- LIHC.CMP[cancer.samples, as.character(LIHC.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- LIHC.CEA[cancer.samples, as.character(LIHC.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  LIHC.diff.40[index] <- m2-m1
  LIHC.fold.40[index] <- m2/m1
}
cat("\n")

LIHC.diff.fold <- data.frame(diff.25 = LIHC.diff.25,
                           diff.33 = LIHC.diff.33,
                           diff.40 = LIHC.diff.40,
                           fold.25 = LIHC.fold.25,
                           fold.33 = LIHC.fold.33,
                           fold.40 = LIHC.fold.40)
rownames(LIHC.diff.fold) <- paste(LIHC.linked.probes.genes$probe, LIHC.linked.probes.genes$genes, sep=".")
LIHC.diff.fold[is.nan(as.matrix(LIHC.diff.fold))] <- NA
save(LIHC.diff.fold, file="../Rdata/LIHC/calc/LIHC-diff-fold.Rdata")
quit(save="no")
