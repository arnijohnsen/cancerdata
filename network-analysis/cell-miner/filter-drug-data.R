# Filter is based on
#  - Finding outlier samples with clustering
#  - Finding drugs with too many unique values

cut.threshold <- 250
min.unique    <- 10
library(WGCNA)

load("../Rdata/cell-miner/data/cell-miner-drug.Rdata")

sample.tree <- flashClust(dist(cell.miner.drug), method="average")
par(cex=0.6)
plot(sample.tree, main="Sample clustering to detect outliers", sub="", xlab="",
     cex.lab=1.5, cex.axis=1.5, cex.main=2)

abline(h=cut.threshold, col="red")

clust <- cutreeStatic(sample.tree, cutHeight=cut.threshold, minSize=10)
table(clust)
keep.samples <- (clust==1)

unique.values <- apply(cell.miner.drug, 2, function(x) {length(unique(x))})
keep.drugs <- (unique.values > 10)

cat("Using", sum(keep.samples), "samples of", dim(cell.miner.drug)[1], "\n")
cat("Using", sum(keep.drugs),   "drugs of",   dim(cell.miner.drug)[2], "\n")
cell.miner.drug.filtered <- cell.miner.drug[keep.samples, keep.drugs]
save(cell.miner.drug.filtered, file="../Rdata/cell-miner/data/cell-miner-drug-filtered.Rdata")
