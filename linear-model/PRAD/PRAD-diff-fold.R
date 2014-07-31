#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/PRAD/info/PRAD-linked-probes-genes.Rdata")
load("../Rdata/PRAD/data/PRAD-CMP.Rdata")
load("../Rdata/PRAD/data/PRAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(PRAD.CMP), rownames(PRAD.CEA))

chunk.size <- 1000
total.size <- dim(PRAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

PRAD.diff.25 <- rep(0, total.size)
PRAD.fold.25 <- rep(0, total.size)
PRAD.diff.33 <- rep(0, total.size)
PRAD.fold.33 <- rep(0, total.size)
PRAD.diff.40 <- rep(0, total.size)
PRAD.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  PRAD.diff.25[index] <- m2-m1
  PRAD.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  PRAD.diff.33[index] <- m2-m1
  PRAD.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- PRAD.CMP[cancer.samples, as.character(PRAD.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- PRAD.CEA[cancer.samples, as.character(PRAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  PRAD.diff.40[index] <- m2-m1
  PRAD.fold.40[index] <- m2/m1
}
cat("\n")

PRAD.diff.fold <- data.frame(diff.25 = PRAD.diff.25,
                           diff.33 = PRAD.diff.33,
                           diff.40 = PRAD.diff.40,
                           fold.25 = PRAD.fold.25,
                           fold.33 = PRAD.fold.33,
                           fold.40 = PRAD.fold.40)
rownames(PRAD.diff.fold) <- paste(PRAD.linked.probes.genes$probe, PRAD.linked.probes.genes$genes, sep=".")
PRAD.diff.fold[is.nan(as.matrix(PRAD.diff.fold))] <- NA
save(PRAD.diff.fold, file="../Rdata/PRAD/calc/PRAD-diff-fold.Rdata")
quit(save="no")
