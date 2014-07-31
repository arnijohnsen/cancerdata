#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/LUAD/info/LUAD-linked-probes-genes.Rdata")
load("../Rdata/LUAD/data/LUAD-CMP.Rdata")
load("../Rdata/LUAD/data/LUAD-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(LUAD.CMP), rownames(LUAD.CEA))

chunk.size <- 1000
total.size <- dim(LUAD.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

LUAD.diff.25 <- rep(0, total.size)
LUAD.fold.25 <- rep(0, total.size)
LUAD.diff.33 <- rep(0, total.size)
LUAD.fold.33 <- rep(0, total.size)
LUAD.diff.40 <- rep(0, total.size)
LUAD.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  LUAD.diff.25[index] <- m2-m1
  LUAD.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  LUAD.diff.33[index] <- m2-m1
  LUAD.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- LUAD.CMP[cancer.samples, as.character(LUAD.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- LUAD.CEA[cancer.samples, as.character(LUAD.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  LUAD.diff.40[index] <- m2-m1
  LUAD.fold.40[index] <- m2/m1
}
cat("\n")

LUAD.diff.fold <- data.frame(diff.25 = LUAD.diff.25,
                           diff.33 = LUAD.diff.33,
                           diff.40 = LUAD.diff.40,
                           fold.25 = LUAD.fold.25,
                           fold.33 = LUAD.fold.33,
                           fold.40 = LUAD.fold.40)
rownames(LUAD.diff.fold) <- paste(LUAD.linked.probes.genes$probe, LUAD.linked.probes.genes$genes, sep=".")
LUAD.diff.fold[is.nan(as.matrix(LUAD.diff.fold))] <- NA
save(LUAD.diff.fold, file="../Rdata/LUAD/calc/LUAD-diff-fold.Rdata")
quit(save="no")
