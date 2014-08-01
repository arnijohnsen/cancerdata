library(WGCNA)
allowWGCNAThreads()
load("../Rdata/cell-miner/data/cell-miner-expression.Rdata")

n <- 10000

start.time <- proc.time()
expression.modules <- blockwiseModules(cell.miner.expression[,1:n],
                        power=6, TOMType="unsigned", minModuleSize=30,
                        reassignThreshold=0, mergeCutHeight=0.25,
                        numericLabels=TRUE, pamRespectsDendro=FALSE,
                        saveTOMs=FALSE,
                        verbose = 3)
print(proc.time()-start.time)
merged.colors <- labels2colors(expression.modules$colors)
plotDendroAndColors(expression.modules$dendrograms[[1]], 
                    merged.colors[expression.modules$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

