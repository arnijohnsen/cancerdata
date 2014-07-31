# Load files
load("../Rdata/LIHC/calc/LIHC-linear-ME.Rdata")
load("../Rdata/LIHC/calc/LIHC-nm-expr.Rdata")
load("../Rdata/LIHC/calc/LIHC-quantiles.Rdata")
load("../Rdata/LIHC/calc/LIHC-diff-fold.Rdata")

LIHC.statistics <- cbind(LIHC.linear.ME, LIHC.nm.expr, LIHC.quantiles, LIHC.diff.fold)
colnames(LIHC.statistics)[6] <- "nm.expr"
save(LIHC.statistics, file="../Rdata/LIHC/calc/LIHC-statistics.Rdata")
quit(save="no")
