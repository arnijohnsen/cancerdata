# Load data files
cat("Loading data files..\n")
load("../Rdata/BRCA/info/BRCA-linked-probes-genes.Rdata")
load("../Rdata/BRCA/data/BRCA-CMP.Rdata")

# Compute diffs
cat("Computing diff and fold 25..\n")
start.time <- proc.time()
diff.fold.25 <- apply(BRCA.CMP, 2, function(x) {
  m1 <- median(x[x>0.25],na.rm=T)
  m2 <- median(x[x<=0.25],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 33..\n")
diff.fold.33 <- apply(BRCA.CMP, 2, function(x) {
  m1 <- median(x[x>0.33],na.rm=T)
  m2 <- median(x[x<=0.33],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 40..\n")
diff.fold.40 <- apply(BRCA.CMP, 2, function(x) {
  m1 <- median(x[x>0.40],na.rm=T)
  m2 <- median(x[x<=0.40],na.rm=T)
  c(m1-m2, m1/m2)
})

BRCA.diff.fold <- data.frame(diff.25 = diff.fold.25[1,as.character(BRCA.linked.probes.genes$probes)],
                             diff.33 = diff.fold.33[1,as.character(BRCA.linked.probes.genes$probes)],
                             diff.40 = diff.fold.40[1,as.character(BRCA.linked.probes.genes$probes)],
                             fold.25 = diff.fold.25[2,as.character(BRCA.linked.probes.genes$probes)],
                             fold.33 = diff.fold.33[2,as.character(BRCA.linked.probes.genes$probes)],
                             fold.40 = diff.fold.40[2,as.character(BRCA.linked.probes.genes$probes)])

rownames(BRCA.diff.fold)<- paste(BRCA.linked.probes.genes$probe,
                                 BRCA.linked.probes.genes$genes, sep=".")

save(BRCA.diff.fold, file="../Rdata/BRCA/calc/BRCA-diff-fold.Rdata")
quit(save="no")
