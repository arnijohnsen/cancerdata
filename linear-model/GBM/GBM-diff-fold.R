#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/GBM/info/GBM-linked-probes-genes.Rdata")
load("../Rdata/GBM/data/GBM-CMP.Rdata")
load("../Rdata/GBM/data/GBM-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(GBM.CMP), rownames(GBM.CEA))

chunk.size <- 1000
total.size <- dim(GBM.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

GBM.diff.25 <- rep(0, total.size)
GBM.fold.25 <- rep(0, total.size)
GBM.diff.33 <- rep(0, total.size)
GBM.fold.33 <- rep(0, total.size)
GBM.diff.40 <- rep(0, total.size)
GBM.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  GBM.diff.25[index] <- m2-m1
  GBM.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  GBM.diff.33[index] <- m2-m1
  GBM.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- GBM.CMP[cancer.samples, as.character(GBM.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- GBM.CEA[cancer.samples, as.character(GBM.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  GBM.diff.40[index] <- m2-m1
  GBM.fold.40[index] <- m2/m1
}
cat("\n")

GBM.diff.fold <- data.frame(diff.25 = GBM.diff.25,
                           diff.33 = GBM.diff.33,
                           diff.40 = GBM.diff.40,
                           fold.25 = GBM.fold.25,
                           fold.33 = GBM.fold.33,
                           fold.40 = GBM.fold.40)
rownames(GBM.diff.fold) <- paste(GBM.linked.probes.genes$probe, GBM.linked.probes.genes$genes, sep=".")
GBM.diff.fold[is.nan(as.matrix(GBM.diff.fold))] <- NA
save(GBM.diff.fold, file="../Rdata/GBM/calc/GBM-diff-fold.Rdata")
quit(save="no")
