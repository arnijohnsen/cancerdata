# Load data files
cat("Loading data files..\n")
load("../Rdata/KIRC/info/KIRC-linked-probes-genes.Rdata")
load("../Rdata/KIRC/data/KIRC-CMP.Rdata")

# Compute diffs
cat("Computing diff and fold 25..\n")
start.time <- proc.time()
diff.fold.25 <- apply(KIRC.CMP, 2, function(x) {
  m1 <- median(x[x>0.25],na.rm=T)
  m2 <- median(x[x<=0.25],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 33..\n")
diff.fold.33 <- apply(KIRC.CMP, 2, function(x) {
  m1 <- median(x[x>0.33],na.rm=T)
  m2 <- median(x[x<=0.33],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 40..\n")
diff.fold.40 <- apply(KIRC.CMP, 2, function(x) {
  m1 <- median(x[x>0.40],na.rm=T)
  m2 <- median(x[x<=0.40],na.rm=T)
  c(m1-m2, m1/m2)
})

KIRC.diff.fold <- data.frame(diff.25 = diff.fold.25[1,as.character(KIRC.linked.probes.genes$probes)],
                             diff.33 = diff.fold.33[1,as.character(KIRC.linked.probes.genes$probes)],
                             diff.40 = diff.fold.40[1,as.character(KIRC.linked.probes.genes$probes)],
                             fold.25 = diff.fold.25[2,as.character(KIRC.linked.probes.genes$probes)],
                             fold.33 = diff.fold.33[2,as.character(KIRC.linked.probes.genes$probes)],
                             fold.40 = diff.fold.40[2,as.character(KIRC.linked.probes.genes$probes)])

rownames(KIRC.diff.fold)<- paste(KIRC.linked.probes.genes$probe,
                                 KIRC.linked.probes.genes$genes, sep=".")

save(KIRC.diff.fold, file="../Rdata/KIRC/calc/KIRC-diff-fold.Rdata")
quit(save="no")
