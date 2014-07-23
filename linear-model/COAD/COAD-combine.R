# Load files
load("../Rdata/COAD/calc/COAD-linear-ME.Rdata")
load("../Rdata/COAD/calc/COAD-nm-expr.Rdata")
load("../Rdata/COAD/calc/COAD-quantiles.Rdata")

COAD.statistics <- cbind(COAD.linear.ME, COAD.nm.expr, COAD.quantiles)
colnames(COAD.statistics)[6] <- "nm.expr"
save(COAD.statistics, file="../Rdata/COAD/calc/COAD-statistics.Rdata")
quit(save="no")
