library(hgu133plus2.db)

load("../Rdata/cell-miner/info/cell-miner-affy-probes.Rdata")

x <- hgu133plus2SYMBOL

cell.miner.gene.list <- unlist(as.list(x[cell.miner.affy.probes]))

save(cell.miner.gene.list, file="../Rdata/cell-miner/info/cell-miner-gene-list.Rdata")
quit(save="no")
