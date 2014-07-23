# Load files
load("../Rdata/LIHC/calc/LIHC-linear-ME.Rdata")
load("../Rdata/LIHC/calc/LIHC-nm-expr.Rdata")
load("../Rdata/LIHC/calc/LIHC-quantiles.Rdata")

LIHC.statistics <- cbind(LIHC.linear.ME, LIHC.nm.expr, LIHC.quantiles)
colnames(LIHC.statistics)[6] <- "nm.expr"
save(LIHC.statistics, file="../Rdata/LIHC/calc/LIHC-statistics.Rdata")
quit(save="no")
