library(WGCNA)

cell.miner.raw <- read.csv("../rawdata/cell-miner/drug/DTP_NCI60.csv", 
                            header           = TRUE,
                            sep              = ";",
                            row.names        = 2,
                            stringsAsFactors = FALSE)

rownames(cell.miner.raw) <- paste("NSC", rownames(cell.miner.raw), sep="")

cell.miner.unfiltered <- t(cell.miner.raw[ ,-c(1:4, 65:66)])

gsg <- goodSamplesGenes(cell.miner.unfiltered, verbose=3)

cell.miner.drug    <- cell.miner.unfiltered[gsg$goodSamples, gsg$goodGenes]
cell.miner.drug.list     <- colnames(cell.miner.drug)
cell.miner.drug.samples  <- rownames(cell.miner.drug)
cell.miner.drug.metadata <- cell.miner.raw[cell.miner.drug.list ,c(1:4)]

save(cell.miner.drug,          file="../Rdata/cell-miner/data/cell-miner-drug.Rdata")
save(cell.miner.drug.metadata, file="../Rdata/cell-miner/info/cell-miner-drug-metadata.Rdata")
save(cell.miner.drug.list,     file="../Rdata/cell-miner/info/cell-miner-drug-list.Rdata")
save(cell.miner.drug.samples,  file="../Rdata/cell-miner/info/cell-miner-drug-samples.Rdata")

quit(save="no")
