#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")
load("../Rdata/BRCA/data/BRCA-CMP.Rdata")
load("../Rdata/BRCA/data/BRCA-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(BRCA.CMP), rownames(BRCA.CEA))

chunk.size <- 1000
total.size <- dim(BRCA.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

BRCA.diff.25 <- rep(0, total.size)
BRCA.fold.25 <- rep(0, total.size)
BRCA.diff.33 <- rep(0, total.size)
BRCA.fold.33 <- rep(0, total.size)
BRCA.diff.40 <- rep(0, total.size)
BRCA.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  BRCA.diff.25[index] <- m2-m1
  BRCA.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  BRCA.diff.33[index] <- m2-m1
  BRCA.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- BRCA.CMP[cancer.samples, as.character(BRCA.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- BRCA.CEA[cancer.samples, as.character(BRCA.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  BRCA.diff.40[index] <- m2-m1
  BRCA.fold.40[index] <- m2/m1
}
cat("\n")

BRCA.diff.fold <- data.frame(diff.25 = BRCA.diff.25,
                           diff.33 = BRCA.diff.33,
                           diff.40 = BRCA.diff.40,
                           fold.25 = BRCA.fold.25,
                           fold.33 = BRCA.fold.33,
                           fold.40 = BRCA.fold.40)
rownames(BRCA.diff.fold) <- paste(BRCA.linked.probes.genes$probe, BRCA.linked.probes.genes$genes, sep=".")
BRCA.diff.fold[is.nan(as.matrix(BRCA.diff.fold))] <- NA
save(BRCA.diff.fold, file="../Rdata/BRCA/calc/BRCA-diff-fold.Rdata")
quit(save="no")
