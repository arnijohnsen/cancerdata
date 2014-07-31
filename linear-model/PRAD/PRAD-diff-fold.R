# Load data files
cat("Loading data files..\n")
load("../Rdata/PRAD/info/PRAD-linked-probes-genes.Rdata")
load("../Rdata/PRAD/data/PRAD-CMP.Rdata")

# Compute diffs
cat("Computing diff and fold 25..\n")
start.time <- proc.time()
diff.fold.25 <- apply(PRAD.CMP, 2, function(x) {
  m1 <- median(x[x>0.25],na.rm=T)
  m2 <- median(x[x<=0.25],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 33..\n")
diff.fold.33 <- apply(PRAD.CMP, 2, function(x) {
  m1 <- median(x[x>0.33],na.rm=T)
  m2 <- median(x[x<=0.33],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 40..\n")
diff.fold.40 <- apply(PRAD.CMP, 2, function(x) {
  m1 <- median(x[x>0.40],na.rm=T)
  m2 <- median(x[x<=0.40],na.rm=T)
  c(m1-m2, m1/m2)
})

PRAD.diff.fold <- data.frame(diff.25 = diff.fold.25[1,as.character(PRAD.linked.probes.genes$probes)],
                             diff.33 = diff.fold.33[1,as.character(PRAD.linked.probes.genes$probes)],
                             diff.40 = diff.fold.40[1,as.character(PRAD.linked.probes.genes$probes)],
                             fold.25 = diff.fold.25[2,as.character(PRAD.linked.probes.genes$probes)],
                             fold.33 = diff.fold.33[2,as.character(PRAD.linked.probes.genes$probes)],
                             fold.40 = diff.fold.40[2,as.character(PRAD.linked.probes.genes$probes)])

rownames(PRAD.diff.fold)<- paste(PRAD.linked.probes.genes$probe,
                                 PRAD.linked.probes.genes$genes, sep=".")

save(PRAD.diff.fold, file="../Rdata/PRAD/calc/PRAD-diff-fold.Rdata")
quit(save="no")
