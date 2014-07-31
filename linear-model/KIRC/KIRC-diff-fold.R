#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/KIRC/info/KIRC-linked-probes-genes.Rdata")
load("../Rdata/KIRC/data/KIRC-CMP.Rdata")
load("../Rdata/KIRC/data/KIRC-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(KIRC.CMP), rownames(KIRC.CEA))

chunk.size <- 1000
total.size <- dim(KIRC.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

KIRC.diff.25 <- rep(0, total.size)
KIRC.fold.25 <- rep(0, total.size)
KIRC.diff.33 <- rep(0, total.size)
KIRC.fold.33 <- rep(0, total.size)
KIRC.diff.40 <- rep(0, total.size)
KIRC.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  KIRC.diff.25[index] <- m2-m1
  KIRC.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  KIRC.diff.33[index] <- m2-m1
  KIRC.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- KIRC.CMP[cancer.samples, as.character(KIRC.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- KIRC.CEA[cancer.samples, as.character(KIRC.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  KIRC.diff.40[index] <- m2-m1
  KIRC.fold.40[index] <- m2/m1
}
cat("\n")

KIRC.diff.fold <- data.frame(diff.25 = KIRC.diff.25,
                           diff.33 = KIRC.diff.33,
                           diff.40 = KIRC.diff.40,
                           fold.25 = KIRC.fold.25,
                           fold.33 = KIRC.fold.33,
                           fold.40 = KIRC.fold.40)
rownames(KIRC.diff.fold) <- paste(KIRC.linked.probes.genes$probe, KIRC.linked.probes.genes$genes, sep=".")
KIRC.diff.fold[is.nan(as.matrix(KIRC.diff.fold))] <- NA
save(KIRC.diff.fold, file="../Rdata/KIRC/calc/KIRC-diff-fold.Rdata")
quit(save="no")
