# Load files
load("../Rdata/READ/calc/READ-linear-ME.Rdata")
load("../Rdata/READ/calc/READ-nm-expr.Rdata")
load("../Rdata/READ/calc/READ-quantiles.Rdata")
load("../Rdata/READ/calc/READ-diff-fold.Rdata")

READ.statistics <- cbind(READ.linear.ME, READ.nm.expr, READ.quantiles, READ.diff.fold)
colnames(READ.statistics)[6] <- "nm.expr"
save(READ.statistics, file="../Rdata/READ/calc/READ-statistics.Rdata")
quit(save="no")
