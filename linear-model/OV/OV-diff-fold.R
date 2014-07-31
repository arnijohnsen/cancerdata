#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/OV/info/OV-linked-probes-genes.Rdata")
load("../Rdata/OV/data/OV-CMP.Rdata")
load("../Rdata/OV/data/OV-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(OV.CMP), rownames(OV.CEA))

chunk.size <- 1000
total.size <- dim(OV.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

OV.diff.25 <- rep(0, total.size)
OV.fold.25 <- rep(0, total.size)
OV.diff.33 <- rep(0, total.size)
OV.fold.33 <- rep(0, total.size)
OV.diff.40 <- rep(0, total.size)
OV.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  OV.diff.25[index] <- m2-m1
  OV.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  OV.diff.33[index] <- m2-m1
  OV.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- OV.CMP[cancer.samples, as.character(OV.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- OV.CEA[cancer.samples, as.character(OV.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  OV.diff.40[index] <- m2-m1
  OV.fold.40[index] <- m2/m1
}
cat("\n")

OV.diff.fold <- data.frame(diff.25 = OV.diff.25,
                           diff.33 = OV.diff.33,
                           diff.40 = OV.diff.40,
                           fold.25 = OV.fold.25,
                           fold.33 = OV.fold.33,
                           fold.40 = OV.fold.40)
rownames(OV.diff.fold) <- paste(OV.linked.probes.genes$probe, OV.linked.probes.genes$genes, sep=".")
OV.diff.fold[is.nan(as.matrix(OV.diff.fold))] <- NA
save(OV.diff.fold, file="../Rdata/OV/calc/OV-diff-fold.Rdata")
quit(save="no")
