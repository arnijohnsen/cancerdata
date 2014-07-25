library(WGCNA)

cell.miner.raw <- read.csv("../rawdata/cell-miner/expression/NCI60_RNA__Affy_HG_U133_Plus_2.0_GCRMA.txt",
                            header           = TRUE,
                            sep              = "\t",
                            row.names        = 2,
                            stringsAsFactors = FALSE)

cell.miner.unfiltered         <- t(cell.miner.raw[ ,-c(1:3)])

gsg <- goodSamplesGenes(cell.miner.unfiltered, verbose=3)

cell.miner.expression         <- cell.miner.unfiltered[gsg$goodSamples, gsg$goodGenes]
cell.miner.affy.probes        <- colnames(cell.miner.expression)
cell.miner.expression.samples <- rownames(cell.miner.expression)

save(cell.miner.expression,          file="../Rdata/cell-miner/data/cell-miner-expression.Rdata")
save(cell.miner.affy.probes,         file="../Rdata/cell-miner/info/cell-miner-affy-probes.Rdata")
save(cell.miner.expression.samples,  file="../Rdata/cell-miner/info/cell-miner-expression-samples.Rdata")

quit(save="no")
