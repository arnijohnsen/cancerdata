#Load library
library(matrixStats)

# Load data files
cat("Loading data files\n")
load("../Rdata/READ/info/READ-linked-probes-genes.Rdata")
load("../Rdata/READ/data/READ-CMP.Rdata")
load("../Rdata/READ/data/READ-CEA.Rdata")


# Resize matrices and use only intersecting samples
cancer.samples <- intersect(rownames(READ.CMP), rownames(READ.CEA))

chunk.size <- 1000
total.size <- dim(READ.linked.probes.genes)[1]
chunks     <- ceiling(total.size/chunk.size)

READ.diff.25 <- rep(0, total.size)
READ.fold.25 <- rep(0, total.size)
READ.diff.33 <- rep(0, total.size)
READ.fold.33 <- rep(0, total.size)
READ.diff.40 <- rep(0, total.size)
READ.fold.40 <- rep(0, total.size)

cat("Running loop for 25\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] > 0.25
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] <= 0.25
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  READ.diff.25[index] <- m2-m1
  READ.fold.25[index] <- m2/m1
}
cat("\n")

cat("Running loop for 33\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] > 0.33
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] <= 0.33
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  READ.diff.33[index] <- m2-m1
  READ.fold.33[index] <- m2/m1
}
cat("\n")
cat("Running loop for 40\n")
pb <- txtProgressBar(min=1, max=chunks, style=3)
for(i in 1:chunks){
  setTxtProgressBar(pb, i)
  index <- (1+chunk.size*(i-1)):(min(chunk.size*i, total.size))
  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] > 0.40
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  m1 <-  colMedians(as.matrix(y), na.rm=T)
  x <- READ.CMP[cancer.samples, as.character(READ.linked.probes.genes$probes[index])] <= 0.40
  x[x==0] <- NA
  y <- READ.CEA[cancer.samples, as.character(READ.linked.probes.genes$genes[index])]*x
  m2 <-  colMedians(as.matrix(y), na.rm=T)
  READ.diff.40[index] <- m2-m1
  READ.fold.40[index] <- m2/m1
}
cat("\n")

READ.diff.fold <- data.frame(diff.25 = READ.diff.25,
                           diff.33 = READ.diff.33,
                           diff.40 = READ.diff.40,
                           fold.25 = READ.fold.25,
                           fold.33 = READ.fold.33,
                           fold.40 = READ.fold.40)
rownames(READ.diff.fold) <- paste(READ.linked.probes.genes$probe, READ.linked.probes.genes$genes, sep=".")
READ.diff.fold[is.nan(as.matrix(READ.diff.fold))] <- NA
save(READ.diff.fold, file="../Rdata/READ/calc/READ-diff-fold.Rdata")
quit(save="no")
