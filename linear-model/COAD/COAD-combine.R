# Load files
load("../Rdata/COAD/calc/COAD-linear-ME.Rdata")
load("../Rdata/COAD/calc/COAD-nm-expr.Rdata")
load("../Rdata/COAD/calc/COAD-quantiles.Rdata")
load("../Rdata/COAD/calc/COAD-diff-fold.Rdata")

COAD.statistics <- cbind(COAD.linear.ME, COAD.nm.expr, COAD.quantiles, COAD.diff.fold)
colnames(COAD.statistics)[6] <- "nm.expr"
save(COAD.statistics, file="../Rdata/COAD/calc/COAD-statistics.Rdata")
quit(save="no")
