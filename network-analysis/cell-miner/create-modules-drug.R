library(WGCNA)
allowWGCNAThreads()
load("../Rdata/cell-miner/data/cell-miner-drug-filtered.Rdata")

n <- 10000

start.time <- proc.time()
drug.modules <- blockwiseModules(cell.miner.drug.filtered[,1:n],
                        power=6, TOMType="unsigned", minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=TRUE, pamRespectsDendro=FALSE,
                        saveTOMs=FALSE, saveTOMFileBase="drugExpressionTOM",
                        verbose = 3)
print(proc.time()-start.time)
merged.colors <- labels2colors(drug.modules$colors)
plotDendroAndColors(drug.modules$dendrograms[[1]], merged.colors[drug.modules$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

