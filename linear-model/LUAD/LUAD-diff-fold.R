# Load data files
cat("Loading data files..\n")
load("../Rdata/LUAD/info/LUAD-linked-probes-genes.Rdata")
load("../Rdata/LUAD/data/LUAD-CMP.Rdata")

# Compute diffs
cat("Computing diff and fold 25..\n")
start.time <- proc.time()
diff.fold.25 <- apply(LUAD.CMP, 2, function(x) {
  m1 <- median(x[x>0.25],na.rm=T)
  m2 <- median(x[x<=0.25],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 33..\n")
diff.fold.33 <- apply(LUAD.CMP, 2, function(x) {
  m1 <- median(x[x>0.33],na.rm=T)
  m2 <- median(x[x<=0.33],na.rm=T)
  c(m1-m2, m1/m2)
})
cat("Computing diff and fold 40..\n")
diff.fold.40 <- apply(LUAD.CMP, 2, function(x) {
  m1 <- median(x[x>0.40],na.rm=T)
  m2 <- median(x[x<=0.40],na.rm=T)
  c(m1-m2, m1/m2)
})

LUAD.diff.fold <- data.frame(diff.25 = diff.fold.25[1,as.character(LUAD.linked.probes.genes$probes)],
                             diff.33 = diff.fold.33[1,as.character(LUAD.linked.probes.genes$probes)],
                             diff.40 = diff.fold.40[1,as.character(LUAD.linked.probes.genes$probes)],
                             fold.25 = diff.fold.25[2,as.character(LUAD.linked.probes.genes$probes)],
                             fold.33 = diff.fold.33[2,as.character(LUAD.linked.probes.genes$probes)],
                             fold.40 = diff.fold.40[2,as.character(LUAD.linked.probes.genes$probes)])

rownames(LUAD.diff.fold)<- paste(LUAD.linked.probes.genes$probe,
                                 LUAD.linked.probes.genes$genes, sep=".")

save(LUAD.diff.fold, file="../Rdata/LUAD/calc/LUAD-diff-fold.Rdata")
quit(save="no")
